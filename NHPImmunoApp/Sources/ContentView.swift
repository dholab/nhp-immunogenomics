import SwiftUI

enum SidebarItem: String, Hashable, CaseIterable {
    case alleles = "Alleles"
    case provisional = "Provisional"
    case submissions = "Submissions"

    var icon: String {
        switch self {
        case .alleles: "list.bullet.rectangle"
        case .provisional: "flask"
        case .submissions: "arrow.triangle.pull"
        }
    }
}

struct ContentView: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(AlleleViewModel.self) private var alleleVM
    @Environment(ProvisionalViewModel.self) private var provisionalVM
    @State private var selection: SidebarItem? = .alleles

    var body: some View {
        NavigationSplitView {
            List(SidebarItem.allCases, id: \.self, selection: $selection) { item in
                Label(item.rawValue, systemImage: item.icon)
            }
            .navigationTitle("NHP Immuno")
        } detail: {
            if settings.isConfigured {
                switch selection {
                case .alleles:
                    AlleleListView()
                case .provisional:
                    ProvisionalListView()
                case .submissions:
                    SubmissionListView()
                case nil:
                    ContentUnavailableView("Select a section",
                        systemImage: "sidebar.left",
                        description: Text("Choose a section from the sidebar."))
                }
            } else {
                ContentUnavailableView {
                    Label("Setup Required", systemImage: "gear")
                } description: {
                    Text("Open Settings (\u{2318},) to configure your GitHub token.")
                }
            }
        }
        .task {
            guard settings.isConfigured, alleleVM.alleles.isEmpty else { return }
            await alleleVM.load(owner: settings.owner, repo: settings.repo)
        }
        .onChange(of: settings.connectionStatus) { _, newStatus in
            if case .connected = newStatus {
                Task { await alleleVM.load(owner: settings.owner, repo: settings.repo) }
            }
        }
    }
}
