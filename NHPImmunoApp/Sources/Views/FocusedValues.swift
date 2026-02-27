import SwiftUI

/// Focused value keys for menu-to-view communication.
extension FocusedValues {
    /// An async closure that refreshes data for the currently active tab.
    @Entry var refreshAction: (() async -> Void)? = nil
}
