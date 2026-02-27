import SwiftUI

enum SidebarItem: String, Hashable, CaseIterable {
    case provisional = "Provisional"
    case submissions = "Submissions"

    var icon: String {
        switch self {
        case .provisional: "flask"
        case .submissions: "arrow.triangle.pull"
        }
    }
}

struct ContentView: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var provisionalVM
    @State private var selection: SidebarItem? = .provisional

    var body: some View {
        NavigationSplitView {
            List(SidebarItem.allCases, id: \.self, selection: $selection) { item in
                Label(item.rawValue, systemImage: item.icon)
            }
            .navigationTitle("NHP Immunogenomics")
        } detail: {
            if settings.isConfigured {
                switch selection {
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
    }
}
