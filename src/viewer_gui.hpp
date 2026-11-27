#pragma once

#include <GLFW/glfw3.h>

struct ViewerState;
class Geometry;

/// Manages Dear ImGui context lifecycle, HUD overlay rendering, and user parameter bindings for molecular viewing.
class ViewerGui {
public:
    /// Constructs an uninitialized viewer GUI manager.
    ViewerGui() = default;

    /// Releases ImGui backend and context resources.
    ~ViewerGui();

    /// Initializes Dear ImGui context and attaches GLFW window callbacks and OpenGL3 pipeline.
    bool init(GLFWwindow* window);

    /// Detaches Dear ImGui render backend and releases memory allocations.
    void shutdown();

    /// Initiates a new ImGui input polling and frame rendering cycle.
    void newFrame();

    /// Constructs interactive control widgets, trajectory timelines, and rendering configuration panels.
    void render(ViewerState& state);

    /// Dispatches generated ImGui geometry draw commands to the active OpenGL rasterizer.
    void draw();

    /// Inverts the graphical user interface visibility flag.
    void toggleVisible() { show_controls = !show_controls; }

    /// Toggles the modal file selection dialog for opening a molecular structure or trajectory file.
    void toggleFileDialog();

    /// Opens the modal file selection dialog for choosing a new molecular trajectory file.
    void openFileDialog();

    /// Closes the modal file selection dialog if currently displayed.
    void closeFileDialog();

    /// Checks whether the modal file selection dialog is currently active and visible.
    bool isFileDialogOpen() const;

    /// Toggles the modal file selection dialog for choosing the output destination of a PNG screenshot.
    void toggleScreenshotDialog();

    /// Opens the modal file selection dialog for choosing the output destination of a PNG screenshot.
    void openScreenshotDialog();

    /// Closes the modal screenshot destination file selection dialog if currently displayed.
    void closeScreenshotDialog();

    /// Checks whether the modal screenshot save dialog is currently active and visible.
    bool isScreenshotDialogOpen() const;

    /// Opens modal file selection dialog to choose output destination for current molecular geometry in XYZ format.
    void openSaveGeometryDialog();

    /// Closes modal geometry destination file selection dialog if currently displayed.
    void closeSaveGeometryDialog();

    /// Checks whether modal geometry destination file selection dialog is currently active.
    bool isSaveGeometryDialogOpen() const;

    /// Opens modal file selection dialog to choose output destination for full molecular trajectory in XYZ format.
    void openSaveTrajectoryDialog();

    /// Closes modal trajectory destination file selection dialog if currently displayed.
    void closeSaveTrajectoryDialog();

    /// Checks whether modal trajectory destination file selection dialog is currently active.
    bool isSaveTrajectoryDialogOpen() const;

    /// Returns whether any part of the graphical user interface overlay is currently displayed.
    bool isVisible(const ViewerState* state = nullptr) const;

    /// Queries whether the Dear ImGui subsystem currently consumes cursor pointing events.
    bool wantsCaptureMouse() const;

    /// Queries whether the Dear ImGui subsystem currently consumes keyboard text input.
    bool wantsCaptureKeyboard() const;

    /// Queries whether an ImGui text input field is actively capturing keyboard characters.
    bool wantsTextInput() const;

private:
    bool initialized   = false;
    bool show_controls = false;

    /// Constructs parameter tuning and camera optics control panels.
    void renderControls(ViewerState& state);

    /// Manages structure loading and screenshot export file browser modals.
    void renderFileDialogs(ViewerState& state);

    /// Constructs the geometric measurement and selected atom inspector overlay.
    void renderMeasurements(ViewerState& state);

    /// Renders in-viewport floating atom index and symbol identification badges.
    void renderFloatingLabels(const ViewerState& state, const Geometry& geom);
};
