#include "viewer_gui.hpp"
#include "viewer_state.hpp"

#include <ImGuiFileDialog.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <cmath>
#include <string>

ViewerGui::~ViewerGui() {
    shutdown();
}

bool ViewerGui::init(GLFWwindow* window) {
    IMGUI_CHECKVERSION(); ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;

    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true) ) return false;
    if (!ImGui_ImplOpenGL3_Init("#version 330 core")) return false;

    initialized = true;

    return true;
}

void ViewerGui::shutdown() {
    if (!initialized) return;

    ImGui_ImplOpenGL3_Shutdown(); ImGui_ImplGlfw_Shutdown(); ImGui::DestroyContext();

    initialized = false;
}

void ViewerGui::newFrame() {
    ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame();
}

void ViewerGui::toggleFileDialog() {
    if (isFileDialogOpen()) {
        closeFileDialog();
    }

    else {
        openFileDialog();
    }
}

void ViewerGui::openFileDialog() {
    IGFD::FileDialogConfig config;

    config.path = "."; config.countSelectionMax = 1;

    ImGuiFileDialog::Instance()->OpenDialog("xyz", "Open Molecular Structure / Trajectory", ".xyz,.*", config);
}

void ViewerGui::closeFileDialog() {
    if (isFileDialogOpen()) {
        ImGuiFileDialog::Instance()->Close();
    }
}

bool ViewerGui::isFileDialogOpen() const {
    if (!initialized) return false;
    return ImGuiFileDialog::Instance()->IsOpened("xyz");
}

void ViewerGui::openScreenshotDialog() {
    IGFD::FileDialogConfig config;

    config.path = "."; config.fileName = "screenshot.png"; config.countSelectionMax = 1;

    config.flags = ImGuiFileDialogFlags_ConfirmOverwrite;

    ImGuiFileDialog::Instance()->OpenDialog("screenshot", "Save Molecular View Screenshot", ".png", config);
}

void ViewerGui::closeScreenshotDialog() {
    if (isScreenshotDialogOpen()) {
        ImGuiFileDialog::Instance()->Close();
    }
}

bool ViewerGui::isScreenshotDialogOpen() const {
    if (!initialized) return false;

    return ImGuiFileDialog::Instance()->IsOpened("screenshot");
}

bool ViewerGui::isVisible() const {
    return show_controls || isFileDialogOpen() || isScreenshotDialogOpen();
}

void ViewerGui::render(ViewerState& state) {
    if (!isVisible()) return;

    if (show_controls) {
        ImGui::SetNextWindowSize(ImVec2(360, 480), ImGuiCond_FirstUseEver);

        if (ImGui::Begin("Molecular Viewer Controls", &show_controls)) {
            int total = state.trajectory.frameCount();
            int cur = state.trajectory.currentFrame();

            if (ImGui::CollapsingHeader("Trajectory Playback", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (total > 1) {
                    int slider_frame = cur + 1; float wait_ms = state.trajectory.getWaitMs();

                    if (ImGui::SliderInt("Frame", &slider_frame, 1, total)) {
                        state.trajectory.setFrame(slider_frame - 1); state.trajectory.setPaused(true);
                    }

                    if (ImGui::Button(state.trajectory.isPaused() ? "  Play  " : " Pause ")) {
                        state.trajectory.togglePause();
                    }

                    ImGui::SameLine();

                    if (ImGui::Button(" < Prev ")) {
                        state.trajectory.stepFrame(-1); state.trajectory.setPaused(true);
                    }

                    ImGui::SameLine();

                    if (ImGui::Button(" Next > ")) {
                        state.trajectory.stepFrame(1); state.trajectory.setPaused(true);
                    }

                    int fps = static_cast<int>(std::round(1000 / wait_ms));

                    if (ImGui::SliderInt("Speed (FPS)", &fps, 1, 120)) {
                        state.trajectory.adjustSpeed((1000 / fps) / wait_ms);
                    }
                }

                else {
                    ImGui::TextDisabled("Static structure");
                }
            }

            if (ImGui::CollapsingHeader("Display Style", ImGuiTreeNodeFlags_DefaultOpen)) {
                const char* styles[] = { "Ball and Stick", "Space Filling", "Wireframe" };

                int style_idx = static_cast<int>(state.style);

                if (ImGui::Combo("Style", &style_idx, styles, 3)) {
                    state.style = static_cast<ViewStyle>(style_idx);
                }

                ImGui::Checkbox("Dark Background", &state.dark_mode);

                if (ImGui::Button("Save Screenshot (S)")) {
                    openScreenshotDialog();
                }
            }

            if (ImGui::CollapsingHeader("Camera Orbit", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (ImGui::Button("Reset Camera (R)")) {
                    state.camera.reset();
                }

                ImGui::Text("Distance: %.2f", state.camera.distance);

                ImGui::Text("Target: (%.2f, %.2f, %.2f)", state.camera.target.x, state.camera.target.y, state.camera.target.z);
            }
        }

        ImGui::End();
    }

    if (isFileDialogOpen()) {
        ImVec2 min_size = ImVec2(500,  350);
        ImVec2 max_size = ImVec2(1000, 700);

        if (ImGuiFileDialog::Instance()->Display("xyz", ImGuiWindowFlags_NoCollapse, min_size, max_size)) {
            if (ImGuiFileDialog::Instance()->IsOk()) {
                state.loadTrajectory(ImGuiFileDialog::Instance()->GetFilePathName());
            }

            ImGuiFileDialog::Instance()->Close();
        }
    }

    if (isScreenshotDialogOpen()) {
        ImVec2 min_size = ImVec2(500,  350);
        ImVec2 max_size = ImVec2(1000, 700);

        if (ImGuiFileDialog::Instance()->Display("screenshot", ImGuiWindowFlags_NoCollapse, min_size, max_size)) {
            if (ImGuiFileDialog::Instance()->IsOk()) {
                state.pending_screenshot_path = ImGuiFileDialog::Instance()->GetFilePathName();
            }

            ImGuiFileDialog::Instance()->Close();
        }
    }
}

void ViewerGui::draw() {
    ImGui::Render(); ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool ViewerGui::wantsCaptureMouse() const {
    if (!initialized || !isVisible()) return false;

    return ImGui::GetIO().WantCaptureMouse;
}

bool ViewerGui::wantsCaptureKeyboard() const {
    if (!initialized || !isVisible()) return false;

    return ImGui::GetIO().WantCaptureKeyboard;
}

bool ViewerGui::wantsTextInput() const {
    if (!initialized || !isVisible()) return false;

    return ImGui::GetIO().WantTextInput;
}
