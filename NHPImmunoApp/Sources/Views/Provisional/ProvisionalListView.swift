import SwiftUI

struct ProvisionalListView: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var vm
    @State private var showDeleteConfirmation = false

    var body: some View {
        @Bindable var vm = vm
        mainContent()
            .navigationTitle("Provisional Alleles")
            .toolbar { toolbarContent }
            .sheet(isPresented: $vm.showAddSheet) {
                AddProvisionalSheet()
            }
            .sheet(item: editBinding) { allele in
                EditProvisionalSheet(allele: allele)
            }
            .confirmationDialog("Delete Provisional Allele",
                                isPresented: $showDeleteConfirmation,
                                presenting: vm.selectedProvisional) { allele in
                Button("Delete \(allele.name)", role: .destructive) {
                    Task { await vm.deleteAlleles(api: settings.api, names: [allele.name]) }
                }
            } message: { allele in
                Text("This will create a PR to remove \(allele.name). The change requires review before it takes effect.")
            }
            .overlay { loadingOverlay }
            .alert("Error", isPresented: errorBinding) {
                Button("OK") { vm.errorMessage = nil }
            } message: {
                Text(vm.errorMessage ?? "")
            }
            .alert("Success", isPresented: successBinding) {
                Button("OK") { vm.successMessage = nil }
            } message: {
                Text(vm.successMessage ?? "")
            }
            .task {
                await vm.loadProvisionals(api: settings.api)
            }
    }

    /// Binding that presents the edit sheet using the selected allele as the item.
    private var editBinding: Binding<ProvisionalAllele?> {
        .init(
            get: { vm.showEditSheet ? vm.selectedProvisional : nil },
            set: { if $0 == nil { vm.showEditSheet = false } }
        )
    }

    @ViewBuilder
    private func mainContent() -> some View {
        @Bindable var vm = vm
        if vm.provisionals.isEmpty && !vm.isLoading {
            ContentUnavailableView {
                Label("No Provisional Alleles", systemImage: "flask")
            } description: {
                Text("Add provisional allele sequences that are not yet in IPD.")
            } actions: {
                Button("Add Alleles") { vm.showAddSheet = true }
            }
        } else {
            Table(vm.provisionals, selection: $vm.selectedProvisionalID) {
                TableColumn("Name", value: \.name).width(min: 150, ideal: 200)
                TableColumn("Species", value: \.species).width(min: 50, ideal: 60)
                TableColumn("Locus", value: \.locus).width(min: 50, ideal: 60)
                TableColumn("Class") { a in Text(a.alleleClass.isEmpty ? "-" : a.alleleClass) }
                    .width(min: 40, ideal: 50)
                TableColumn("Submitter", value: \.submitter).width(min: 80, ideal: 120)
                TableColumn("Date Added", value: \.dateAdded).width(min: 80, ideal: 100)
                TableColumn("Notes", value: \.notes).width(min: 80, ideal: 150)
            }
        }
    }

    @ToolbarContentBuilder
    private var toolbarContent: some ToolbarContent {
        ToolbarItemGroup {
            Button("Add", systemImage: "plus") { vm.showAddSheet = true }
                .help("Add provisional alleles")
                .keyboardShortcut("n", modifiers: .command)
            Button("Edit", systemImage: "pencil") { vm.showEditSheet = true }
                .disabled(vm.selectedProvisional == nil)
                .help("Edit selected allele")
            Button("Delete", systemImage: "trash") { showDeleteConfirmation = true }
                .disabled(vm.selectedProvisional == nil)
                .help("Delete selected allele")
            Divider()
            Button("Refresh", systemImage: "arrow.clockwise") {
                Task { await vm.loadProvisionals(api: settings.api) }
            }
            .help("Refresh provisional alleles")
            .keyboardShortcut("r", modifiers: .command)
        }
    }

    @ViewBuilder
    private var loadingOverlay: some View {
        if vm.isSubmitting {
            ZStack {
                Color.black.opacity(0.2)
                ProgressView("Submitting...")
                    .padding()
                    .background(.regularMaterial, in: .rect(cornerRadius: 12))
            }
        } else if vm.isLoading {
            ProgressView("Loading...")
        }
    }

    private var errorBinding: Binding<Bool> {
        .init(get: { vm.errorMessage != nil }, set: { if !$0 { vm.errorMessage = nil } })
    }

    private var successBinding: Binding<Bool> {
        .init(get: { vm.successMessage != nil }, set: { if !$0 { vm.successMessage = nil } })
    }
}
