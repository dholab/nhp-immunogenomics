import Foundation

/// Fetch and cache allele data from GitHub Pages.
actor AlleleService {
    private let pagesBase: URL
    private let session = URLSession.shared
    private var cachedIndex: AlleleIndex?

    init(owner: String, repo: String) {
        var components = URLComponents()
        components.scheme = "https"
        components.host = "\(owner).github.io"
        components.path = "/\(repo)/"
        self.pagesBase = components.url ?? URL(string: "https://\(owner).github.io/\(repo)/")!
    }

    /// Fetch alleles.json from GitHub Pages.
    func fetchIndex(forceRefresh: Bool = false) async throws -> AlleleIndex {
        if !forceRefresh, let cached = cachedIndex { return cached }
        let url = pagesBase.appendingPathComponent("alleles.json")
        let (data, response) = try await session.data(from: url)
        guard let http = response as? HTTPURLResponse, http.statusCode == 200 else {
            throw AlleleServiceError.fetchFailed("Failed to fetch alleles.json")
        }
        let index = try JSONDecoder().decode(AlleleIndex.self, from: data)
        cachedIndex = index
        return index
    }

    /// Fetch a GenBank file from GitHub Pages. Returns the raw text.
    func fetchGenBank(project: String, species: String, locus: String) async throws -> String {
        let path = "data/\(project)/\(species)/\(locus).gb"
        let url = pagesBase.appendingPathComponent(path)
        let (data, response) = try await session.data(from: url)
        guard let http = response as? HTTPURLResponse, http.statusCode == 200,
              let text = String(data: data, encoding: .utf8) else {
            throw AlleleServiceError.fetchFailed("Failed to fetch \(path)")
        }
        return text
    }

    func clearCache() {
        cachedIndex = nil
    }
}

enum AlleleServiceError: LocalizedError {
    case fetchFailed(String)
    var errorDescription: String? {
        switch self { case .fetchFailed(let msg): msg }
    }
}
