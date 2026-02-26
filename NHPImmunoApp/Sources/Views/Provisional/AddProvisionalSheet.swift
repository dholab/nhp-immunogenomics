import SwiftUI

struct AddProvisionalSheet: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var vm
    @Environment(\.dismiss) private var dismiss

    @State private var species = ""
    @State private var locus = ""
    @State private var alleleClass = "I"
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
        !species.isEmpty && !locus.isEmpty && !submitter.isEmpty && !fastaContent.isEmpty && validationErrors.isEmpty
    }

    var body: some View {
        NavigationStack {
            Form {
                Section("Allele Metadata") {
                    TextField("Species prefix", text: $species, prompt: Text("e.g., Mamu"))
                    TextField("Locus", text: $locus, prompt: Text("e.g., A1, B, E"))
                    Picker("MHC Class", selection: $alleleClass) {
                        Text("Class I").tag("I")
                        Text("Class II").tag("II")
                        Text("None").tag("")
                    }
                    Picker("Sequence Type", selection: $seqType) {
                        Text("Coding").tag("coding")
                        Text("Genomic").tag("genomic")
                    }
                    TextField("Submitter", text: $submitter, prompt: Text("e.g., J. Karl"))
                    TextField("Notes (optional)", text: $notes)
                }

                Section("FASTA Sequences") {
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

                    if !fastaContent.isEmpty {
                        Text(fastaPreview)
                            .font(.body.monospaced())
                            .foregroundStyle(.secondary)
                            .lineLimit(6)
                            .frame(maxWidth: .infinity, alignment: .leading)
                            .textSelection(.enabled)
                    }

                    if !validationErrors.isEmpty {
                        ForEach(validationErrors, id: \.self) { error in
                            Label(error, systemImage: "exclamationmark.triangle.fill")
                                .foregroundStyle(.red)
                                .font(.caption)
                        }
                    }
                }

                if let error = submissionError {
                    Section {
                        Label(error, systemImage: "xmark.circle.fill")
                            .foregroundStyle(.red)
                    }
                }
            }
            .formStyle(.grouped)
            .navigationTitle("Add Provisional Alleles")
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

    private var fastaPreview: String {
        let lines = fastaContent.split(separator: "\n", maxSplits: 5).map(String.init)
        return lines.prefix(5).joined(separator: "\n") + (lines.count > 5 ? "\n..." : "")
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
                api: settings.api, species: species, locus: locus,
                alleleClass: alleleClass, seqType: seqType,
                submitter: submitter, notes: notes,
                fastaContent: fastaContent
            )
            if let error = vm.errorMessage {
                submissionError = error
            }
        }
    }
}
