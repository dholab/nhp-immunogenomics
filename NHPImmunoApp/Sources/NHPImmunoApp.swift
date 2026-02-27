import SwiftUI

@main
struct NHPImmunoApp: App {
    @State private var settings = SettingsViewModel()
    @State private var alleleVM = AlleleViewModel()
    @State private var provisionalVM = ProvisionalViewModel()

    var body: some Scene {
        WindowGroup {
            ContentView()
                .environment(settings)
                .environment(alleleVM)
                .environment(provisionalVM)
        }
        .defaultSize(width: 1100, height: 700)
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
        ViewMenuCommands()
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

// MARK: - View Menu

private struct ViewMenuCommands: Commands {
    @FocusedValue(\.toggleInspector) private var toggleInspector

    var body: some Commands {
        CommandGroup(after: .toolbar) {
            Button("Toggle Inspector") {
                toggleInspector?()
            }
            .keyboardShortcut("i", modifiers: [.command, .option])
            .disabled(toggleInspector == nil)
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
