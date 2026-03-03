import SwiftUI

enum ProvisionalFilter {
    case mhc, kir
}

struct ProvisionalListView: View {
    let filter: ProvisionalFilter

    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var vm
    @State private var showDeleteConfirmation = false

    private var filteredProvisionals: [ProvisionalAllele] {
        switch filter {
        case .mhc: vm.provisionals.filter { !$0.locus.hasPrefix("KIR") }
        case .kir: vm.provisionals.filter { $0.locus.hasPrefix("KIR") }
        }
    }

    private var navigationTitle: String {
        switch filter {
        case .mhc: "MHC Provisional Alleles"
        case .kir: "KIR Provisional Alleles"
        }
    }

    var body: some View {
        @Bindable var vm = vm
        mainContent()
            .navigationTitle(navigationTitle)
            .toolbar { toolbarContent }
            .sheet(isPresented: $vm.showAddSheet) { addSheet }
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
            .focusedValue(\.refreshAction) {
                await vm.loadProvisionals(api: settings.api)
            }
    }

    @ViewBuilder
    private var addSheet: some View {
        if filter == .kir {
            AddKIRProvisionalSheet()
        } else {
            AddProvisionalSheet()
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
        let emptyLabel = filter == .kir ? "No KIR Provisional Alleles" : "No MHC Provisional Alleles"
        let emptyDescription = filter == .kir
            ? "Add provisional KIR allele sequences that are not yet in IPD."
            : "Add provisional MHC allele sequences that are not yet in IPD."
        if filteredProvisionals.isEmpty && !vm.isLoading {
            ContentUnavailableView {
                Label(emptyLabel, systemImage: "flask")
            } description: {
                Text(emptyDescription)
            } actions: {
                Button("Add Alleles") { vm.showAddSheet = true }
            }
        } else {
            Table(filteredProvisionals, selection: $vm.selectedProvisionalID) {
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

// MARK: - Add KIR Provisional Sheet

struct AddKIRProvisionalSheet: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var vm
    @Environment(\.dismiss) private var dismiss

    @State private var species = ""
    @State private var seqType = "coding"
    @State private var submitter = ""
    @State private var notes = ""
    @State private var fastaContent = ""
    @State private var fastaFileName = ""
    @State private var showFileImporter = false
    @State private var validationErrors: [String] = []
    @State private var submissionError: String?
    @State private var parsedRecordCount = 0

    private var isValid: Bool {
        !species.isEmpty && !submitter.isEmpty && !fastaContent.isEmpty && validationErrors.isEmpty
    }

    var body: some View {
        NavigationStack {
            Form {
                Section("Allele Metadata") {
                    TextField("Species prefix", text: $species, prompt: Text("e.g., Mamu"))
                    Picker("Sequence Type", selection: $seqType) {
                        Text("Coding").tag("coding")
                        Text("Genomic").tag("genomic")
                    }
                    TextField("Submitter", text: $submitter, prompt: Text("e.g., J. Karl"))
                    TextField("Notes (optional)", text: $notes)
                }

                Section {
                    Text("The locus will be auto-detected by searching against all KIR references for the species.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }

                Section("FASTA Sequences") {
                    fastaPickerRow
                    fastaPreviewRow
                    fastaValidationRows
                }

                if let error = submissionError {
                    Section {
                        Label(error, systemImage: "xmark.circle.fill")
                            .foregroundStyle(.red)
                    }
                }
            }
            .formStyle(.grouped)
            .navigationTitle("Add KIR Provisional Alleles")
            .frame(minWidth: 500, minHeight: 400)
            .toolbar {
                ToolbarItem(placement: .cancellationAction) {
                    Button("Cancel") { dismiss() }
                }
                ToolbarItem(placement: .confirmationAction) {
                    Button("Submit") { submit() }
                        .disabled(!isValid || vm.isSubmitting)
                }
            }
            .fileImporter(isPresented: $showFileImporter,
                          allowedContentTypes: [.plainText, .data],
                          allowsMultipleSelection: false) { result in
                handleFileImport(result)
            }
        }
    }

    @ViewBuilder
    private var fastaPickerRow: some View {
        HStack {
            Button("Choose File...") { showFileImporter = true }
            if !fastaFileName.isEmpty {
                Text(fastaFileName)
                    .foregroundStyle(.secondary)
                    .lineLimit(1)
            }
            Spacer()
            if parsedRecordCount > 0 {
                Text("\(parsedRecordCount) sequence(s)")
                    .foregroundStyle(.secondary)
            }
        }
    }

    @ViewBuilder
    private var fastaPreviewRow: some View {
        if !fastaContent.isEmpty {
            let lines = fastaContent.split(separator: "\n", maxSplits: 5).map(String.init)
            let preview = lines.prefix(5).joined(separator: "\n") + (lines.count > 5 ? "\n..." : "")
            Text(preview)
                .font(.body.monospaced())
                .foregroundStyle(.secondary)
                .lineLimit(6)
                .frame(maxWidth: .infinity, alignment: .leading)
                .textSelection(.enabled)
        }
    }

    @ViewBuilder
    private var fastaValidationRows: some View {
        ForEach(validationErrors, id: \.self) { error in
            Label(error, systemImage: "exclamationmark.triangle.fill")
                .foregroundStyle(.red)
                .font(.caption)
        }
    }

    private func handleFileImport(_ result: Result<[URL], Error>) {
        switch result {
        case .success(let urls):
            guard let url = urls.first else { return }
            guard url.startAccessingSecurityScopedResource() else { return }
            defer { url.stopAccessingSecurityScopedResource() }
            do {
                fastaContent = try String(contentsOf: url, encoding: .utf8)
                fastaFileName = url.lastPathComponent
                validate()
            } catch {
                validationErrors = ["Failed to read file: \(error.localizedDescription)"]
            }
        case .failure(let error):
            validationErrors = ["File picker error: \(error.localizedDescription)"]
        }
    }

    private func validate() {
        let records = FASTAParser.parse(fastaContent)
        parsedRecordCount = records.count
        validationErrors = FASTAParser.validateRecords(records)
    }

    private func submit() {
        submissionError = nil
        Task {
            await vm.submitNewAlleles(
                api: settings.api, species: species, locus: "",
                alleleClass: "", seqType: seqType,
                submitter: submitter, notes: notes,
                fastaContent: fastaContent
            )
            if let error = vm.errorMessage {
                submissionError = error
            }
        }
    }
}
