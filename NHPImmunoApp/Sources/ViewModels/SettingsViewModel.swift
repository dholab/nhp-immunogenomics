import Foundation
import Observation

@MainActor
@Observable
final class SettingsViewModel {
    var token: String = ""
    var owner: String = "dholab"
    var repo: String = "nhp-immunogenomics"
    var connectionStatus: ConnectionStatus = .untested
    var isTesting = false

    enum ConnectionStatus: Equatable {
        case untested
        case testing
        case connected(String)
        case failed(String)

        var label: String {
            switch self {
            case .untested: "Not tested"
            case .testing: "Testing..."
            case .connected(let name): "Connected to \(name)"
            case .failed(let msg): "Failed: \(msg)"
            }
        }
    }

    init() {
        token = KeychainService.loadToken() ?? ""
        owner = UserDefaults.standard.string(forKey: "githubOwner") ?? "dholab"
        repo = UserDefaults.standard.string(forKey: "githubRepo") ?? "nhp-immunogenomics"
    }

    /// Save settings to Keychain and UserDefaults.
    /// Call this explicitly (e.g., on window close or test connection), not on every keystroke.
    func save() {
        do {
            if token.isEmpty {
                KeychainService.deleteToken()
            } else {
                try KeychainService.saveToken(token)
            }
        } catch {
            connectionStatus = .failed(error.localizedDescription)
        }
        UserDefaults.standard.set(owner, forKey: "githubOwner")
        UserDefaults.standard.set(repo, forKey: "githubRepo")
    }

    func testConnection() async {
        save()
        connectionStatus = .testing
        isTesting = true
        defer { isTesting = false }

        do {
            let name = try await api.verifyAccess()
            connectionStatus = .connected(name)
        } catch {
            connectionStatus = .failed(error.localizedDescription)
        }
    }

    var isConfigured: Bool {
        !token.isEmpty && !owner.isEmpty && !repo.isEmpty
    }

    /// Convenience accessor for creating a GitHubAPI instance from current settings.
    var api: GitHubAPI {
        GitHubAPI(owner: owner, repo: repo, token: token)
    }
}
