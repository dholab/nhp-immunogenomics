import SwiftUI

@main
struct NHPImmunoApp: App {
    @State private var settings = SettingsViewModel()
    @State private var provisionalVM = ProvisionalViewModel()

    var body: some Scene {
        WindowGroup {
            ContentView()
                .environment(settings)
                .environment(provisionalVM)
        }
        .defaultSize(width: 900, height: 600)
        .commands {
            appMenus
        }

        Settings {
            SettingsView()
                .environment(settings)
        }
    }

    @CommandsBuilder
    private var appMenus: some Commands {
        FileMenuCommands()
        HelpMenuCommands()
    }
}

// MARK: - File Menu

private struct FileMenuCommands: Commands {
    @FocusedValue(\.refreshAction) private var refreshAction

    var body: some Commands {
        CommandGroup(after: .newItem) {
            Button("Refresh Data") {
                guard let action = refreshAction else { return }
                Task { await action() }
            }
            .keyboardShortcut("r", modifiers: .command)
            .disabled(refreshAction == nil)
        }
    }
}

// MARK: - Help Menu

private struct HelpMenuCommands: Commands {
    var body: some Commands {
        CommandGroup(replacing: .help) {
            Link("NHP Immunogenomics on GitHub",
                 destination: URL(string: "https://github.com/dholab/nhp-immunogenomics")!)
        }
    }
}
