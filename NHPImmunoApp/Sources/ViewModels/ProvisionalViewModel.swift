import Foundation
import Observation

@MainActor
@Observable
final class ProvisionalViewModel {
    // Data
    var provisionals: [ProvisionalAllele] = []
    var submissions: [PullRequestInfo] = []

    // State
    var isLoading = false
    var isLoadingSubmissions = false
    var isSubmitting = false
    var errorMessage: String?
    var successMessage: String?
    var selectedProvisionalID: ProvisionalAllele.ID?

    /// The currently selected provisional allele, resolved from the ID.
    var selectedProvisional: ProvisionalAllele? {
        guard let id = selectedProvisionalID else { return nil }
        return provisionals.first { $0.id == id }
    }

    // Sheet presentation
    var showAddSheet = false
    var showEditSheet = false

    func loadProvisionals(api: GitHubAPI) async {
        isLoading = true
        errorMessage = nil
        defer { isLoading = false }

        let service = ProvisionalService(api: api)
        do {
            provisionals = try await service.fetchManifest()
        } catch {
            errorMessage = error.localizedDescription
        }
    }

    func loadSubmissions(api: GitHubAPI) async {
        isLoadingSubmissions = true
        defer { isLoadingSubmissions = false }
        do {
            submissions = try await api.listPullRequests(headPrefix: "provisional/")
            for i in submissions.indices {
                let pr = submissions[i]
                if let head = try? await api.getBranchSHA(pr.headRef) {
                    submissions[i].checksStatus = (try? await api.getCheckStatus(sha: head)) ?? .unknown
                }
            }
        } catch {
            errorMessage = "Failed to load submissions: \(error.localizedDescription)"
        }
    }

    func submitNewAlleles(api: GitHubAPI, species: String, locus: String,
                          alleleClass: String, seqType: String, submitter: String,
                          notes: String, fastaContent: String) async {
        isSubmitting = true
        errorMessage = nil
        successMessage = nil
        defer { isSubmitting = false }

        let service = ProvisionalService(api: api)
        do {
            let prURL = try await service.submitNewAlleles(
                species: species, locus: locus, alleleClass: alleleClass,
                seqType: seqType, submitter: submitter, notes: notes,
                fastaContent: fastaContent
            )
            successMessage = "PR created: \(prURL.absoluteString)"
            showAddSheet = false
            await loadSubmissions(api: api)
        } catch {
            errorMessage = error.localizedDescription
        }
    }

    func editAllele(api: GitHubAPI, allele: ProvisionalAllele) async {
        isSubmitting = true
        errorMessage = nil
        successMessage = nil
        defer { isSubmitting = false }

        let service = ProvisionalService(api: api)
        do {
            let prURL = try await service.editAllele(allele)
            successMessage = "Edit PR created: \(prURL.absoluteString)"
            showEditSheet = false
            await loadSubmissions(api: api)
        } catch {
            errorMessage = error.localizedDescription
        }
    }

    func deleteAlleles(api: GitHubAPI, names: Set<String>) async {
        isSubmitting = true
        errorMessage = nil
        successMessage = nil
        defer { isSubmitting = false }

        let service = ProvisionalService(api: api)
        do {
            let prURL = try await service.deleteAlleles(names)
            successMessage = "Delete PR created: \(prURL.absoluteString)"
            await loadProvisionals(api: api)
            await loadSubmissions(api: api)
        } catch {
            errorMessage = error.localizedDescription
        }
    }
}
