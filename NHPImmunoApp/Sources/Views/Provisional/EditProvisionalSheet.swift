import SwiftUI

struct EditProvisionalSheet: View {
    @Environment(SettingsViewModel.self) private var settings
    @Environment(ProvisionalViewModel.self) private var vm
    @Environment(\.dismiss) private var dismiss

    @State var allele: ProvisionalAllele
    @State private var submissionError: String?

    var body: some View {
        NavigationStack {
            Form {
                Section("Allele") {
                    LabeledContent("Name", value: allele.name)
                    LabeledContent("Species", value: allele.species)
                    LabeledContent("Locus", value: allele.locus)
                    LabeledContent("Class", value: allele.alleleClass.isEmpty ? "-" : allele.alleleClass)
                    LabeledContent("Seq Type", value: allele.seqType)
                    LabeledContent("Date Added", value: allele.dateAdded)
                }

                Section("Editable Fields") {
                    TextField("Submitter", text: $allele.submitter)
                    TextField("Notes", text: $allele.notes, axis: .vertical)
                        .lineLimit(3...6)
                }

                if let error = submissionError {
                    Section {
                        Label(error, systemImage: "xmark.circle.fill")
                            .foregroundStyle(.red)
                    }
                }
            }
            .formStyle(.grouped)
            .navigationTitle("Edit \(allele.name)")
            .frame(minWidth: 450, minHeight: 300)
            .toolbar {
                ToolbarItem(placement: .cancellationAction) {
                    Button("Cancel") { dismiss() }
                }
                ToolbarItem(placement: .confirmationAction) {
                    Button("Save") { save() }
                        .disabled(vm.isSubmitting)
                }
            }
        }
    }

    private func save() {
        submissionError = nil
        Task {
            await vm.editAllele(api: settings.api, allele: allele)
            if let error = vm.errorMessage {
                submissionError = error
            }
        }
    }
}
