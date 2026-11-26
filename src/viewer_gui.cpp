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

    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    io.IniFilename              = nullptr;
    io.ConfigNavCaptureKeyboard =   false;

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

void ViewerGui::toggleScreenshotDialog() {
    if (isScreenshotDialogOpen()) {
        closeScreenshotDialog();
    }

    else {
        openScreenshotDialog();
    }
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
        ImGui::SetNextWindowSize(ImVec2(380, 560), ImGuiCond_FirstUseEver);

        if (ImGui::Begin("Molecular Viewer Controls", &show_controls)) {
            if (ImGui::Button("Open Molecule / Trajectory (O)", ImVec2(-1.0f, 0.0f))) {
                openFileDialog();
            }

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

                    int fps = static_cast<int>(std::round(1000.0f / wait_ms));

                    if (ImGui::SliderInt("Speed (FPS)", &fps, 1, 120)) {
                        state.trajectory.adjustSpeed((1000.0f / static_cast<float>(fps)) / wait_ms);
                    }

                    const char* loop_modes[] = { "Loop", "Ping-Pong", "Once" };

                    int mode_idx = static_cast<int>(state.trajectory.getLoopMode());

                    if (ImGui::Combo("Loop Mode", &mode_idx, loop_modes, 3)) {
                        state.trajectory.setLoopMode(static_cast<LoopMode>(mode_idx));
                    }
                }

                else {
                    ImGui::TextDisabled("Static structure");
                }
            }

            if (ImGui::CollapsingHeader("Display & Meshing", ImGuiTreeNodeFlags_DefaultOpen)) {
                const char* styles[] = { "Ball and Stick", "Space Filling", "Wireframe" };

                int style_idx = static_cast<int>(state.style);

                if (ImGui::Combo("Style", &style_idx, styles, 3)) {
                    state.style = static_cast<ViewStyle>(style_idx);
                }

                if (ImGui::ColorEdit3("Background Color", &state.bg_color.x)) {
                    state.dark_mode = (state.bg_color.x + state.bg_color.y + state.bg_color.z) < 1.5f;
                }

                ImGui::Checkbox("Transparent Background", &state.transparent_bg);

                ImGui::SliderFloat("Atom Scale", &state.atom_scale, 0.2f, 3.0f, "%.2f");

                if (ImGui::Checkbox("Smooth Shading", &state.smooth_shading)) {
                    state.mesh_dirty = true;
                }

                if (ImGui::SliderInt("Sphere Subdiv", &state.sphere_subdivisions, 0, 4)) {
                    state.mesh_dirty = true;
                }

                if (ImGui::SliderInt("Cylinder Sectors", &state.cylinder_sectors, 4, 48)) {
                    state.mesh_dirty = true;
                }

                if (ImGui::Button("Save Screenshot (S)")) {
                    openScreenshotDialog();
                }
            }

            if (ImGui::CollapsingHeader("Lighting & Shading", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::SliderFloat("Ambient",   &state.ambient_intensity,  0.0f,   1.0f, "%.2f");
                ImGui::SliderFloat("Diffuse",   &state.diffuse_intensity,  0.0f,   1.0f, "%.2f");
                ImGui::SliderFloat("Specular",  &state.specular_intensity, 0.0f,   1.0f, "%.2f");
                ImGui::SliderFloat("Shininess", &state.shininess,          1.0f, 128.0f, "%.1f");

                ImGui::SliderFloat("Light Offset X", &state.light_offset.x, -2.0f, 2.0f, "%.2f");
                ImGui::SliderFloat("Light Offset Y", &state.light_offset.y, -2.0f, 2.0f, "%.2f");

                if (ImGui::Button("Reset Lighting")) {
                    state.ambient_intensity  = 0.35f;
                    state.diffuse_intensity  = 0.55f;
                    state.specular_intensity = 0.25f;
                    state.shininess          = 32.0f;

                    state.light_offset = vec3(0.4f, 0.7f, 0.0f);
                }
            }

            if (ImGui::CollapsingHeader("Bonding & Connectivity", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (ImGui::SliderFloat("Binding Factor", &state.binding_factor, 0.005f, 0.035f, "%.4f")) {
                    state.trajectory.rebind(state.binding_factor, state.bond_radius);
                }

                if (ImGui::SliderFloat("Bond Radius", &state.bond_radius, 0.01f, 0.30f, "%.3f")) {
                    state.trajectory.rebind(state.binding_factor, state.bond_radius);
                }

                ImGui::ColorEdit3("Bond Color", &state.bond_color.x);

                if (ImGui::Button("Reset Bonds")) {
                    state.binding_factor = 0.013f;
                    state.bond_radius    = 0.090f;

                    state.bond_color = vec3(0.65f, 0.65f, 0.65f);

                    state.trajectory.rebind(state.binding_factor, state.bond_radius);
                }
            }

            if (ImGui::CollapsingHeader("Camera & Optics", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (ImGui::Button("Reset Camera (R)")) {
                    state.camera.reset();
                }

                const char* proj_modes[] = { "Perspective", "Orthographic" };

                int proj_idx = static_cast<int>(state.camera.projection);

                if (ImGui::Combo("Projection (P)", &proj_idx, proj_modes, 2)) {
                    state.camera.projection = static_cast<ProjectionMode>(proj_idx);
                }

                ImGui::SliderFloat("Field of View", &state.camera.fov_deg, 10.0f,   120.0f, "%.1f deg");
                ImGui::SliderFloat("Near Plane",    &state.camera.z_near,   0.01f,   10.0f, "%.2f"    );
                ImGui::SliderFloat("Far Plane",     &state.camera.z_far,  100.00f, 5000.0f, "%.0f"    );

                ImGui::Text("Distance: %.2f", state.camera.distance);

                ImGui::Text("Target: (%.2f, %.2f, %.2f)", state.camera.target.x, state.camera.target.y, state.camera.target.z);
            }
        }

        ImGui::End();
    }

    if (isFileDialogOpen()) {
        ImVec2 min_size = ImVec2( 500.0f, 350.0f);
        ImVec2 max_size = ImVec2(1000.0f, 700.0f);

        if (ImGuiFileDialog::Instance()->Display("xyz", ImGuiWindowFlags_NoCollapse, min_size, max_size)) {
            if (ImGuiFileDialog::Instance()->IsOk()) {
                state.loadTrajectory(ImGuiFileDialog::Instance()->GetFilePathName());
            }

            ImGuiFileDialog::Instance()->Close();
        }
    }

    if (isScreenshotDialogOpen()) {
        ImVec2 min_size = ImVec2( 500.0f, 350.0f);
        ImVec2 max_size = ImVec2(1000.0f, 700.0f);

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
