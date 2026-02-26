import Foundation

/// Minimal pull request information for the submissions view.
struct PullRequestInfo: Identifiable, Sendable {
    let number: Int
    let title: String
    let state: String
    let url: URL
    let createdAt: Date
    let headRef: String
    var checksStatus: ChecksStatus

    var id: Int { number }

    enum ChecksStatus: Sendable {
        case pending
        case success
        case failure
        case unknown
    }

    /// Pre-formatted relative date string for display.
    var relativeDate: String {
        let formatter = RelativeDateTimeFormatter()
        formatter.unitsStyle = .short
        return formatter.localizedString(for: createdAt, relativeTo: Date())
    }

    /// Parse a PR creation date from ISO 8601 string.
    static func parseDate(_ string: String) -> Date {
        let formatter = ISO8601DateFormatter()
        return formatter.date(from: string) ?? Date()
    }
}
