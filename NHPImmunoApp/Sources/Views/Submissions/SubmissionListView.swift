import SwiftUI

struct SubmissionListView: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var vm

    var body: some View {
        Group {
            if vm.isLoadingSubmissions {
                ProgressView("Loading submissions...")
            } else if vm.submissions.isEmpty {
                ContentUnavailableView {
                    Label("No Open Submissions", systemImage: "checkmark.circle")
                } description: {
                    Text("Provisional allele submissions that are awaiting review will appear here.")
                }
            } else {
                submissionList
            }
        }
        .navigationTitle("Submissions")
        .toolbar {
            ToolbarItem {
                Button("Refresh", systemImage: "arrow.clockwise") {
                    Task { await vm.loadSubmissions(api: settings.api) }
                }
                .help("Refresh submissions")
            }
        }
        .task {
            await vm.loadSubmissions(api: settings.api)
        }
        .focusedValue(\.refreshAction) {
            await vm.loadSubmissions(api: settings.api)
        }
    }

    private var submissionList: some View {
        List(vm.submissions) { pr in
            SubmissionRow(pr: pr)
        }
    }
}

private struct SubmissionRow: View {
    let pr: PullRequestInfo

    var body: some View {
        HStack {
            statusIcon

            VStack(alignment: .leading, spacing: 2) {
                Text(pr.title)
                    .fontWeight(.medium)
                Text("#\(pr.number) opened \(pr.relativeDate)")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }

            Spacer()

            statusLabel

            Link(destination: pr.url) {
                Image(systemName: "arrow.up.right.square")
            }
            .help("Open in GitHub")
        }
        .padding(.vertical, 4)
    }

    @ViewBuilder
    private var statusIcon: some View {
        switch pr.checksStatus {
        case .pending:
            Image(systemName: "clock.fill").foregroundStyle(.orange)
        case .success:
            Image(systemName: "checkmark.circle.fill").foregroundStyle(.green)
        case .failure:
            Image(systemName: "xmark.circle.fill").foregroundStyle(.red)
        case .unknown:
            Image(systemName: "questionmark.circle").foregroundStyle(.secondary)
        }
    }

    @ViewBuilder
    private var statusLabel: some View {
        let (text, color): (String, Color) = switch pr.checksStatus {
        case .pending: ("Pending", .orange)
        case .success: ("Passed", .green)
        case .failure: ("Failed", .red)
        case .unknown: ("Unknown", .secondary)
        }
        Text(text)
            .font(.caption)
            .padding(.horizontal, 6)
            .padding(.vertical, 2)
            .background(color.opacity(0.15), in: .capsule)
            .foregroundStyle(color)
    }
}
