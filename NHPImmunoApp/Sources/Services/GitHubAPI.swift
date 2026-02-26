import Foundation

/// HTTP client for the GitHub REST API.
/// All methods are async and throw ``GitHubAPIError`` on failure.
actor GitHubAPI {
    let owner: String
    let repo: String
    private let token: String
    private let session: URLSession

    init(owner: String, repo: String, token: String) {
        self.owner = owner
        self.repo = repo
        self.token = token
        self.session = URLSession.shared
    }

    // MARK: - URL Building

    /// Build a GitHub API URL from path segments and optional query items.
    /// Each segment is individually percent-encoded to prevent injection.
    private func apiURL(path: String, queryItems: [URLQueryItem] = []) throws -> URL {
        var components = URLComponents()
        components.scheme = "https"
        components.host = "api.github.com"
        components.path = path
        if !queryItems.isEmpty {
            components.queryItems = queryItems
        }
        guard let url = components.url else {
            throw GitHubAPIError.invalidURL(path)
        }
        return url
    }

    /// Build the standard repo API path prefix.
    private var repoPath: String { "/repos/\(owner)/\(repo)" }

    // MARK: - Request Helpers

    private func request(_ path: String, method: String = "GET",
                         body: [String: Any]? = nil,
                         queryItems: [URLQueryItem] = [],
                         accept: String = "application/vnd.github+json") async throws -> (Data, HTTPURLResponse) {
        let url = try apiURL(path: path, queryItems: queryItems)
        var req = URLRequest(url: url)
        req.httpMethod = method
        req.setValue("Bearer \(token)", forHTTPHeaderField: "Authorization")
        req.setValue(accept, forHTTPHeaderField: "Accept")
        req.setValue("2022-11-28", forHTTPHeaderField: "X-GitHub-Api-Version")
        if let body {
            req.httpBody = try JSONSerialization.data(withJSONObject: body)
            req.setValue("application/json", forHTTPHeaderField: "Content-Type")
        }
        let (data, response) = try await session.data(for: req)
        guard let http = response as? HTTPURLResponse else {
            throw GitHubAPIError.invalidResponse
        }
        if !(200...299).contains(http.statusCode) {
            let message = (try? JSONSerialization.jsonObject(with: data) as? [String: Any])?["message"] as? String
            throw GitHubAPIError.httpError(http.statusCode, message ?? String(data: data, encoding: .utf8) ?? "")
        }
        return (data, http)
    }

    private func fetchJSON<T: Decodable>(_ path: String, queryItems: [URLQueryItem] = []) async throws -> T {
        let (data, _) = try await request(path, queryItems: queryItems)
        return try JSONDecoder().decode(T.self, from: data)
    }

    // MARK: - Codable Response Types

    private struct RepoResponse: Decodable {
        let fullName: String
        let permissions: Permissions?

        enum CodingKeys: String, CodingKey {
            case fullName = "full_name"
            case permissions
        }

        struct Permissions: Decodable {
            let push: Bool?
        }
    }

    private struct FileContentResponse: Decodable {
        let content: String?
        let sha: String
    }

    private struct GitRefResponse: Decodable {
        let object: GitObject

        struct GitObject: Decodable {
            let sha: String
        }
    }

    private struct CommitFileResponse: Decodable {
        let commit: CommitInfo

        struct CommitInfo: Decodable {
            let sha: String
        }
    }

    private struct PRResponse: Decodable {
        let number: Int
        let htmlURL: String

        enum CodingKeys: String, CodingKey {
            case number
            case htmlURL = "html_url"
        }
    }

    private struct PRListItem: Decodable {
        let number: Int
        let title: String
        let state: String
        let htmlURL: String
        let createdAt: String
        let head: Head

        enum CodingKeys: String, CodingKey {
            case number, title, state
            case htmlURL = "html_url"
            case createdAt = "created_at"
            case head
        }

        struct Head: Decodable {
            let ref: String
        }
    }

    private struct CheckRunsResponse: Decodable {
        let checkRuns: [CheckRun]

        enum CodingKeys: String, CodingKey {
            case checkRuns = "check_runs"
        }

        struct CheckRun: Decodable {
            let status: String
            let conclusion: String?
        }
    }

    // MARK: - Repository

    /// Verify the token is valid and has access to the repo.
    func verifyAccess() async throws -> String {
        let response: RepoResponse = try await fetchJSON("\(repoPath)")
        return response.fullName
    }

    // MARK: - Contents API

    /// Fetch a file's content and SHA from the repo.
    func getFileContent(path filePath: String, ref: String? = nil) async throws -> (content: String, sha: String) {
        var queryItems: [URLQueryItem] = []
        if let ref { queryItems.append(URLQueryItem(name: "ref", value: ref)) }
        let response: FileContentResponse = try await fetchJSON(
            "\(repoPath)/contents/\(filePath)", queryItems: queryItems)
        guard let encodedContent = response.content else {
            throw GitHubAPIError.unexpectedFormat("Missing content in response")
        }
        let cleaned = encodedContent.replacingOccurrences(of: "\n", with: "")
        guard let decoded = Data(base64Encoded: cleaned),
              let text = String(data: decoded, encoding: .utf8) else {
            throw GitHubAPIError.unexpectedFormat("Unable to decode base64 content")
        }
        return (text, response.sha)
    }

    /// Check if a file exists. Returns its SHA if found, nil otherwise.
    func getFileSHA(path: String, ref: String? = nil) async throws -> String? {
        do {
            let (_, sha) = try await getFileContent(path: path, ref: ref)
            return sha
        } catch GitHubAPIError.httpError(404, _) {
            return nil
        }
    }

    // MARK: - Git Refs (Branches)

    /// Get the SHA of a branch head.
    func getBranchSHA(_ branch: String) async throws -> String {
        let response: GitRefResponse = try await fetchJSON(
            "\(repoPath)/git/ref/heads/\(branch)")
        return response.object.sha
    }

    /// Create a new branch from a SHA.
    func createBranch(name: String, fromSHA sha: String) async throws {
        _ = try await request("\(repoPath)/git/refs", method: "POST", body: [
            "ref": "refs/heads/\(name)",
            "sha": sha,
        ])
    }

    // MARK: - Commits (Contents API)

    /// Create or update a file on a branch. Returns the new commit SHA.
    @discardableResult
    func commitFile(path filePath: String, content: Data, message: String,
                    branch: String, existingSHA: String? = nil) async throws -> String {
        var body: [String: Any] = [
            "message": message,
            "content": content.base64EncodedString(),
            "branch": branch,
        ]
        if let sha = existingSHA {
            body["sha"] = sha
        }
        let (data, _) = try await request("\(repoPath)/contents/\(filePath)",
                                           method: "PUT", body: body)
        let response = try JSONDecoder().decode(CommitFileResponse.self, from: data)
        return response.commit.sha
    }

    /// Delete a file on a branch.
    func deleteFile(path filePath: String, message: String, branch: String, sha: String) async throws {
        _ = try await request("\(repoPath)/contents/\(filePath)", method: "DELETE", body: [
            "message": message,
            "sha": sha,
            "branch": branch,
        ])
    }

    // MARK: - Pull Requests

    /// Create a pull request. Returns the PR number and URL.
    func createPullRequest(title: String, body: String, head: String,
                           base: String = "main", draft: Bool = false) async throws -> (number: Int, url: URL) {
        let (data, _) = try await request("\(repoPath)/pulls", method: "POST", body: [
            "title": title,
            "body": body,
            "head": head,
            "base": base,
            "draft": draft,
        ])
        let response = try JSONDecoder().decode(PRResponse.self, from: data)
        guard let url = URL(string: response.htmlURL) else {
            throw GitHubAPIError.unexpectedFormat("Invalid PR URL")
        }
        return (response.number, url)
    }

    /// List open pull requests, optionally filtered by head branch prefix.
    func listPullRequests(headPrefix: String? = nil) async throws -> [PullRequestInfo] {
        let queryItems = [
            URLQueryItem(name: "state", value: "open"),
            URLQueryItem(name: "per_page", value: "30"),
            URLQueryItem(name: "sort", value: "created"),
            URLQueryItem(name: "direction", value: "desc"),
        ]
        let (data, _) = try await request("\(repoPath)/pulls", queryItems: queryItems)
        let items = try JSONDecoder().decode([PRListItem].self, from: data)
        return items.compactMap { pr -> PullRequestInfo? in
            guard let url = URL(string: pr.htmlURL) else { return nil }
            if let prefix = headPrefix, !pr.head.ref.hasPrefix(prefix) { return nil }
            return PullRequestInfo(
                number: pr.number, title: pr.title, state: pr.state,
                url: url, createdAt: PullRequestInfo.parseDate(pr.createdAt),
                headRef: pr.head.ref, checksStatus: .unknown)
        }
    }

    /// Get check run status for a commit SHA.
    func getCheckStatus(sha: String) async throws -> PullRequestInfo.ChecksStatus {
        let response: CheckRunsResponse = try await fetchJSON(
            "\(repoPath)/commits/\(sha)/check-runs")
        let runs = response.checkRuns
        if runs.isEmpty { return .unknown }
        let statuses = runs.map(\.status)
        let conclusions = runs.compactMap(\.conclusion)
        if statuses.contains("in_progress") || statuses.contains("queued") { return .pending }
        if conclusions.contains("failure") || conclusions.contains("cancelled") { return .failure }
        if conclusions.allSatisfy({ $0 == "success" || $0 == "skipped" }) { return .success }
        return .pending
    }

    // MARK: - Workflow Dispatch

    /// Trigger a workflow_dispatch event on a branch.
    func triggerWorkflow(workflowFileName: String, ref: String,
                         inputs: [String: String] = [:]) async throws {
        var body: [String: Any] = ["ref": ref]
        if !inputs.isEmpty { body["inputs"] = inputs }
        _ = try await request(
            "\(repoPath)/actions/workflows/\(workflowFileName)/dispatches",
            method: "POST", body: body)
    }
}

// MARK: - Errors

enum GitHubAPIError: LocalizedError {
    case invalidResponse
    case invalidURL(String)
    case httpError(Int, String)
    case unexpectedFormat(String)

    var errorDescription: String? {
        switch self {
        case .invalidResponse: "Invalid HTTP response"
        case .invalidURL(let path): "Invalid API URL: \(path)"
        case .httpError(let code, let msg): "GitHub API error \(code): \(msg)"
        case .unexpectedFormat(let msg): msg
        }
    }
}
