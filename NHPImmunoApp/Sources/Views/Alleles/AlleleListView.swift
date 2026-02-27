import SwiftUI

struct AlleleListView: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(AlleleViewModel.self) private var vm
    @State private var sortOrder = [KeyPathComparator(\Allele.name)]
    @State private var showInspector = false

    var body: some View {
        @Bindable var vm = vm

        alleleTable()
            .navigationTitle("Allele Database")
            .searchable(text: $vm.searchText, prompt: "Search alleles...")
            .toolbar { toolbarContent }
            .inspector(isPresented: $showInspector) {
                AlleleDetailView()
                    .inspectorColumnWidth(min: 250, ideal: 300, max: 400)
            }
            .overlay { overlayContent }
            .safeAreaInset(edge: .bottom) { statusBar }
            .focusedValue(\.refreshAction) {
                await vm.load(owner: settings.owner, repo: settings.repo)
            }
            .focusedValue(\.toggleInspector) {
                showInspector.toggle()
            }
            .onChange(of: vm.searchText) { _, _ in vm.applyFilters() }
            .onChange(of: vm.selectedSpecies) { _, _ in vm.applyFilters() }
            .onChange(of: vm.selectedLocus) { _, _ in vm.applyFilters() }
            .onChange(of: vm.selectedClass) { _, _ in vm.applyFilters() }
            .onChange(of: vm.selectedProject) { _, _ in vm.applyFilters() }
            .onChange(of: vm.showProvisionalOnly) { _, _ in vm.applyFilters() }
    }

    private func alleleTable() -> some View {
        @Bindable var vm = vm
        let sorted = vm.filteredAlleles.sorted(using: sortOrder)
        return Table(sorted, selection: $vm.selectedAlleleID, sortOrder: $sortOrder) {
            TableColumn("Accession", value: \.accession) { allele in
                Text(allele.accession).font(.body.monospaced())
            }
            .width(min: 80, ideal: 100)

            TableColumn("Name", value: \.name) { allele in
                AlleleName(allele: allele)
            }
            .width(min: 180, ideal: 250)

            TableColumn("Species", value: \.species)
                .width(min: 50, ideal: 60)

            TableColumn("Locus", value: \.locus)
                .width(min: 50, ideal: 70)

            TableColumn("Class", value: \.alleleClass) { allele in
                Text(allele.alleleClass.isEmpty ? "-" : allele.alleleClass)
            }
            .width(min: 40, ideal: 50)

            TableColumn("Project", value: \.project) { allele in
                Text(allele.projectLabel)
            }
            .width(min: 40, ideal: 50)

            TableColumn("Date Modified", value: \.dateModified)
                .width(min: 80, ideal: 100)
        }
    }

    @ToolbarContentBuilder
    private var toolbarContent: some ToolbarContent {
        ToolbarItem(placement: .principal) {
            AlleleFilterBar()
        }
        ToolbarItem {
            Button("Inspector", systemImage: "sidebar.right") {
                showInspector.toggle()
            }
            .help("Toggle inspector panel")
        }
        ToolbarItem {
            Button("Refresh", systemImage: "arrow.clockwise") {
                Task { await vm.load(owner: settings.owner, repo: settings.repo) }
            }
            .help("Refresh allele data")
        }
    }

    @ViewBuilder
    private var overlayContent: some View {
        if vm.isLoading {
            ProgressView("Loading alleles...")
        } else if let error = vm.errorMessage {
            ContentUnavailableView {
                Label("Error", systemImage: "exclamationmark.triangle")
            } description: {
                Text(error)
            } actions: {
                Button("Retry") {
                    Task { await vm.load(owner: settings.owner, repo: settings.repo) }
                }
            }
        } else if vm.filteredAlleles.isEmpty && !vm.alleles.isEmpty {
            ContentUnavailableView.search(text: vm.searchText)
        }
    }

    private var statusBar: some View {
        HStack {
            Text(vm.alleleCount)
                .foregroundStyle(.secondary)
                .font(.caption)
            Spacer()
            if !vm.mhcVersion.isEmpty {
                Text("MHC \(vm.mhcVersion) / KIR \(vm.nhkirVersion)")
                    .foregroundStyle(.secondary)
                    .font(.caption)
            }
        }
        .padding(.horizontal)
        .padding(.vertical, 6)
        .background(.bar)
    }
}

private struct AlleleName: View {
    let allele: Allele
    var body: some View {
        HStack(spacing: 4) {
            Text(allele.name)
            if allele.isProvisional {
                Text("PROV")
                    .font(.caption2.weight(.semibold))
                    .padding(.horizontal, 4)
                    .padding(.vertical, 1)
                    .background(.orange.opacity(0.2), in: .capsule)
                    .foregroundStyle(.orange)
            }
        }
    }
}
