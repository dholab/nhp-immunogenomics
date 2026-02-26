import SwiftUI

struct AlleleFilterBar: View {
    @Environment(AlleleViewModel.self) private var vm

    var body: some View {
        @Bindable var vm = vm

        HStack(spacing: 12) {
            Picker("Project", selection: $vm.selectedProject) {
                Text("All Projects").tag(String?.none)
                Divider()
                Text("MHC").tag(String?.some("MHC"))
                Text("KIR").tag(String?.some("NHKIR"))
            }
            .frame(width: 130)

            Picker("Species", selection: $vm.selectedSpecies) {
                Text("All Species").tag(String?.none)
                Divider()
                ForEach(vm.availableSpecies, id: \.self) { sp in
                    Text(sp).tag(String?.some(sp))
                }
            }
            .frame(width: 130)

            Picker("Locus", selection: $vm.selectedLocus) {
                Text("All Loci").tag(String?.none)
                Divider()
                ForEach(vm.availableLoci, id: \.self) { locus in
                    Text(locus).tag(String?.some(locus))
                }
            }
            .frame(width: 120)

            Picker("Class", selection: $vm.selectedClass) {
                Text("All Classes").tag(String?.none)
                Divider()
                Text("Class I").tag(String?.some("I"))
                Text("Class II").tag(String?.some("II"))
            }
            .frame(width: 120)

            Toggle("Provisional", isOn: $vm.showProvisionalOnly)

            Spacer()

            if vm.hasActiveFilters {
                Button("Clear Filters", systemImage: "xmark.circle") {
                    vm.clearFilters()
                }
                .buttonStyle(.borderless)
            }
        }
    }
}
