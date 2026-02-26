import SwiftUI

struct SettingsView: View {
    @Environment(SettingsViewModel.self) private var vm

    var body: some View {
        @Bindable var vm = vm

        Form {
            Section("GitHub Authentication") {
                SecureField("Personal Access Token", text: $vm.token)
                    .help("Token needs 'repo' scope for read/write access.")
                Text("Create a token at github.com/settings/tokens with **repo** scope.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }

            Section("Repository") {
                TextField("Owner", text: $vm.owner)
                TextField("Repository", text: $vm.repo)
            }

            Section("Connection") {
                HStack {
                    statusIcon
                    Text(vm.connectionStatus.label)
                        .foregroundStyle(statusColor)
                    Spacer()
                    Button("Test Connection") {
                        Task { await vm.testConnection() }
                    }
                    .disabled(vm.token.isEmpty || vm.isTesting)
                }
            }
        }
        .formStyle(.grouped)
        .frame(width: 450)
        .padding()
        .onDisappear { vm.save() }
    }

    @ViewBuilder
    private var statusIcon: some View {
        switch vm.connectionStatus {
        case .untested:
            Image(systemName: "circle")
                .foregroundStyle(.secondary)
        case .testing:
            ProgressView()
                .controlSize(.small)
        case .connected:
            Image(systemName: "checkmark.circle.fill")
                .foregroundStyle(.green)
        case .failed:
            Image(systemName: "xmark.circle.fill")
                .foregroundStyle(.red)
        }
    }

    private var statusColor: Color {
        switch vm.connectionStatus {
        case .connected: .green
        case .failed: .red
        default: .secondary
        }
    }
}
