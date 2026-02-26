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

        Settings {
            SettingsView()
                .environment(settings)
        }
    }
}
