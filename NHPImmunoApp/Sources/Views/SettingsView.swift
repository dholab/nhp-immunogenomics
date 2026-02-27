import SwiftUI

struct SettingsView: View {
    @Environment(SettingsViewModel.self) private var vm

    var body: some View {
        @Bindable var vm = vm

        Form {
            Section("GitHub Authentication") {
                SecureField("Personal Access Token", text: $vm.token)
                    .help("Fine-grained token with Contents, Pull requests, Actions, and Workflows read/write access.")
                Text("""
                    Create a **fine-grained** token at GitHub → Settings → Developer settings → \
                    Personal access tokens → Fine-grained tokens. \
                    Grant it access to the repository with these permissions: \
                    **Contents** (Read and write), **Pull requests** (Read and write), \
                    **Actions** (Read and write), and **Workflows** (Read and write).
                    """)
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
