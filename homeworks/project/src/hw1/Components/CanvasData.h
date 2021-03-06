#pragma once

#include <imgui/imgui.h>
#include <UGM/UGM.h>

struct CanvasData {
	std::vector<Ubpa::pointf2> points;
	std::vector<std::vector<ImVec2>> draw_points;
	double dx;

	Ubpa::valf2 scrolling{ 0.f,0.f };
	bool opt_enable_grid{ true };
	bool opt_enable_context_menu{ true };
	bool adding_line{ false };
};

#include "details/CanvasData_AutoRefl.inl"
