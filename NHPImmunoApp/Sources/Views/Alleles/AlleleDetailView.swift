import SwiftUI

struct AlleleDetailView: View {
    @Environment(AlleleViewModel.self) private var vm

    var body: some View {
        if let allele = vm.selectedAllele {
            Form {
                Section("Identity") {
                    LabeledContent("Accession", value: allele.accession)
                    LabeledContent("Name", value: allele.name)
                    LabeledContent("Species", value: speciesLabel(allele.species))
                    LabeledContent("Locus", value: allele.locus)
                    LabeledContent("Class", value: allele.alleleClass.isEmpty ? "-" : allele.alleleClass)
                    LabeledContent("Project", value: allele.projectLabel)
                }

                Section("Dates") {
                    LabeledContent("Date Assigned", value: allele.dateAssigned)
                    LabeledContent("Date Modified", value: allele.dateModified)
                }

                if let prev = allele.previousNames, !prev.isEmpty {
                    Section("Previous Names") {
                        ForEach(prev, id: \.self) { name in
                            Text(name)
                                .font(.body.monospaced())
                        }
                    }
                }

                if allele.isProvisional {
                    Section("Provisional") {
                        if let sub = allele.submitter {
                            LabeledContent("Submitter", value: sub)
                        }
                        if let st = allele.seqType {
                            LabeledContent("Seq Type", value: st)
                        }
                    }
                }
            }
            .formStyle(.grouped)
            .navigationTitle(allele.name)
        } else {
            ContentUnavailableView("No Selection",
                systemImage: "arrow.left.circle",
                description: Text("Select an allele to view details."))
        }
    }

    private func speciesLabel(_ prefix: String) -> String {
        if let species = vm.speciesMap[prefix] {
            return "\(prefix) (\(species.scientificName))"
        }
        return prefix
    }
}
