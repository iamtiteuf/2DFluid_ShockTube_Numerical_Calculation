
#include "Glew_Initialization.h"
#include "shader.h"
#include "Mesh.h"
#include "Camera.h"
#include "ColorPicking.h"
#include "Sphere.h"
#include "PhysicsHandler.h"
#include "implot.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_stdlib.h"
#include "imgui_internal.h"



bool start_Simulation = false;
bool IsWidget;
int simVersion = 3;

Shaders shader; // main shader program
Shaders pickingshader; // picking color shader program

GenSphere c_Sphere;//generate sphre vertices

Color_Picking colorPicking; 

Glew_Initialization sys;

std::vector<Vertex> cube_vertex =
{
	// positions          // colors           // normals         // texUV
 { glm::vec3(-0.3,-0.3,-0.3), glm::vec3(1.0,0.0,0.5), glm::vec3(-1,-1,-1), glm::vec2(0, 0) },
 { glm::vec3(-0.3,-0.3,0.3), glm::vec3(1.0,0.0,0.5), glm::vec3(-1,-1, 1), glm::vec2(1, 0) },
 { glm::vec3(-0.3,0.3,-0.3),  glm::vec3(1.0,0.0,0.5), glm::vec3(-1, 1,-1), glm::vec2(0, 1) },
 { glm::vec3(-0.3,0.3,0.3),  glm::vec3(1.0,0.0,0.5), glm::vec3(-1, 1, 1), glm::vec2(1, 1) },
 { glm::vec3(0.3,-0.3,-0.3),  glm::vec3(1.0,0.0,0.5), glm::vec3(1,-1,-1),  glm::vec2(1, 0) },
 { glm::vec3(0.3,-0.3,0.3),   glm::vec3(1.0,0.0,0.5), glm::vec3(1,-1, 1),  glm::vec2(0, 0) },
 { glm::vec3(0.3,0.3,-0.3),  glm::vec3(1.0,0.0,0.5) , glm::vec3(1, 1,-1),  glm::vec2(1, 1) },
 { glm::vec3(0.3,0.3,0.3),    glm::vec3(1.0,0.0,0.5), glm::vec3(1, 1, 1),  glm::vec2(0, 1) }
};
std::vector<unsigned int> indices =
{
   0, 1, 2, 1, 2, 3, 4, 5, 6, 5, 6, 7, // front 
   1, 5, 3, 5, 3, 7, 0, 4, 2, 4, 2, 6, // back 
   0, 1, 4, 1, 4, 5, 2, 3, 6, 3, 6, 7  // sides
};


std::vector<mesh*> All_Mesh;//all objects
std::vector<Transform> Predic_Transforms;


glm::mat4 CamMatrix; //camera Matrix
PhysicsHandler MySim;

void ManageWidget();
void Render_Scene();
void Render_Picking_Scene();
int  Get_Object(glm::vec4 ColorID);
bool Almost_Equal(glm::vec4 a, glm::vec4 b);
void InitiatePicking();
void ClearUp();
void CheckIfWidget(ImGuiIO& io)
{
	if (io.WantCaptureMouse || ImGui::IsAnyItemHovered())
	{
		IsWidget = true;
	}
	else
	{
		IsWidget = false;
	}
}

bool leftclicked = false;
int Cube_ID = 0;
int Cube2_ID = 0;

int main()
{
	sys.Init_Glfw();

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();
	ImGui_ImplGlfw_InitForOpenGL(MainWindow, true);
	ImGui_ImplOpenGL3_Init("#version 330");
	ImGui::StyleColorsDark();
	ImGuiIO& io = ImGui::GetIO();

	Camera camera;

	shader.Create_Normal_Shaders();
	pickingshader.Create_Picking_Shaders();

	// create custome framebuffer
	colorPicking.Init(sys.width, sys.height, glm::vec3(0, 0, 0));
	std::vector<Vertex> Sphere_vert;
	std::vector<unsigned int> Sphere_ind;
	c_Sphere.generateSphere(0.2f, 10, 10, glm::vec3(1, 0, 0), Sphere_vert, Sphere_ind);

	

	//create cubes and assign unique ColorIDs
	

	

	


	

	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	glEnable(GL_DEPTH_TEST);

	while (!glfwWindowShouldClose(MainWindow))
	{
		Predic_Transforms.clear();
		for (auto obj : All_Mesh)
		{
			Predic_Transforms.push_back(obj->transform);
		}

		
		CheckIfWidget(io);
		glfwPollEvents();
		CamMatrix = camera.Camera_Matrix();
		camera.Camera_Controls(true);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		
		Render_Scene();

		//InitiatePicking();


		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		//ui implementation
		
		ManageWidget();

		
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		
		glfwSwapBuffers(MainWindow);
		

	}
	ImPlot::DestroyContext();
	ClearUp();
	glfwDestroyWindow(MainWindow);
	glfwTerminate();
	
}


void Render_Scene()
{
	if (labH5)
	{
		if (restart == true)
		{
			switch (simVersion)
			{
			case 1:
				MySim.Simulation_2();
				break;
			case 2:
				MySim.Simulation_3();
				break;
			case 3:
				MySim.Simulation_4();
				break;
			default:
				MySim.Simulation_4();
				break;
			}
			//MySim.Simulation_2();
			//MySim.Simulation_3();
			//MySim.Simulation_4();
			restart = false;
		}
		
		
	}
	else
	{
		MySim.Poisseulle();
	}
	
	
	for (int i = 0; i < All_Mesh.size(); i++)
	{
		if (start_Simulation)
		{
			PhysicsHandler::Step(All_Mesh[i]->transform);
		}
		
	}
	for (auto Line : All_Mesh)
	{
		Line->RenderMesh(shader.shaderprogram, CamMatrix);
		if (Line->Line_Attached)
		{
			if (Line->show_line)
			{
				Line->DrawLine(shader.shaderprogram, CamMatrix);
			}
		}
	}
	
}
void Render_Picking_Scene()
{
	for (int i = 0; i < All_Mesh.size(); i++)
	{
		All_Mesh[i]->RenderMesh(pickingshader.shaderprogram, CamMatrix);
	}
}
void InitiatePicking()
{
	if (glfwGetMouseButton(MainWindow, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && !IsWidget)
	{
		if (!leftclicked)
		{
			Render_Scene();
			Render_Picking_Scene();
			//draw cubes in the custom frame buffer
			colorPicking.EnableWriting();
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

			Render_Picking_Scene();
			colorPicking.DisableWriting();
			//read color values from custom frame buffer
			Cube_ID = Get_Object(colorPicking.Read_Pixels());

			All_Mesh[Cube_ID]->UpdateColor(glm::vec3(0.9f,0.9f,0.9f));

			for (int i = 0; i < All_Mesh.size(); i++)
			{
				if (i != Cube_ID)
				{
					All_Mesh[i]->UpdateColor(All_Mesh[i]->Color);
				}
			}
			leftclicked = true;
		}

	}
	else if (glfwGetMouseButton(MainWindow, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE)
	{
		if (leftclicked)
		{
			leftclicked = false;
		}

	}
}
int Get_Object(glm::vec4 ColorID)
{
	int e = 0;
	for (int i = 0; i < All_Mesh.size(); i++)
	{
		if (Almost_Equal(All_Mesh[i]->ColorID, ColorID))
		{
			e = i;
			break;

		}
	}
	return e;

}
bool Almost_Equal(glm::vec4 a, glm::vec4 b)
{
	double epsilon = 0.0001;
	bool check_diff = true;
	for (int i = 0; i < 4; i++)
	{
		if (std::abs(a[i] - b[i]) > epsilon)
		{
			check_diff = false;
		}
	}
	return check_diff;
}
void ClearUp()
{
	for (auto obj : All_Mesh)
	{
		obj->ClearMesh();
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}
void ManageWidget()
{
	const char* name = "LAb H5";
	if (labH5)
	{
		name = "LAb H5";
	}
	else { name = "Physics mech"; }
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("Lab Variant"))
		{
			if (ImGui::Button(name))
			{
				labH5 = !labH5;
			}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}
	
	ImGui::Begin("values");
	if (ImGui::InputInt("Simulation version", &simVersion))
	{
		if (simVersion > 3)
		{
			simVersion = 3;
		}
		else if(simVersion < 1)
		{
			simVersion = 1;
		}
	}
	ImGui::Checkbox("Restart", &restart);
	ImGui::InputInt("N", &myN);
	ImGui::InputInt("tmax", &tmax);
	ImGui::InputFloat("K :", &somek);
	ImGui::InputFloat("dx :", &mydx,0.1,0.1,"%.5f");
	ImGui::InputFloat("dt :", &mydt, 0.1, 0.1, "%.7f");
	ImGui::InputFloat("a :", &pisAcc, 0.1, 0.1, "%.5f");
	ImGui::InputFloat("cappa :", &cappa, 0.1, 0.1, "%.3f");
	ImGui::SameLine();
	ImGui::InputFloat("Hartmann", &Ha);
	ImGui::InputInt("Range", &range);
	ImGui::InputFloat("Length of Pipe", &l_of_pipe);
	ImGui::InputFloat("viscousity", &viscousity);
	ImGui::InputFloat("flowrate", &flow_rate);
	ImGui::InputFloat("pipe raduis", &pipe_raduis);
	ImGui::InputFloat("pipe area", &area_of_pipe);

	ImGui::End();
	ImGui::Begin("My graphs");
	ImGui::SetWindowPos(ImVec2(0, 20));
	ImGui::BeginChild("velocity Graphs", ImVec2(0, 1500), true, ImGuiWindowFlags_MenuBar);
	

	if (labH5)
	{
		ImGui::SeparatorText("velocity Profile");
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);

		if (ImPlot::BeginPlot("velocity Plot",ImVec2(-1,500))) {


			ImPlot::PlotLine(" v1 vs x Plot :", MySim.x1.data(), MySim.u1.data(), MySim.x1.size());
			ImPlot::PlotLine(" v2 vs x Plot :", MySim.x2.data(), MySim.u2.data(), MySim.x2.size());
			ImPlot::PlotLine(" v3 vs x Plot :", MySim.x3.data(), MySim.u3.data(), MySim.x3.size());
			ImPlot::PlotLine(" v4 vs x Plot :", MySim.x4.data(), MySim.u4.data(), MySim.x4.size());
			ImPlot::PlotLine(" v5 vs x Plot :", MySim.xn.data(), MySim.un.data(), MySim.xn.size());
			ImPlot::PlotLine(" vn vs x Plot :", MySim.XN.data(), MySim.uN.data(), MySim.XN.size());
			

			ImPlot::EndPlot();
		}
		ImGui::SeparatorText("pressure Profile");
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);

		if (ImPlot::BeginPlot("pressure Plot", ImVec2(-1, 500))) {


			ImPlot::PlotLine("p1 vs x Plot :", MySim.x1.data(), MySim.p1.data(), MySim.x1.size());
			ImPlot::PlotLine("p2 vs x Plot :", MySim.x2.data(), MySim.p2.data(), MySim.x2.size());
			ImPlot::PlotLine("p3 vs x Plot :", MySim.x3.data(), MySim.p3.data(), MySim.x3.size());
		    ImPlot::PlotLine("p4 vs x Plot :", MySim.x4.data(), MySim.p4.data(), MySim.x4.size());
	    	ImPlot::PlotLine("p4 vs x Plot :", MySim.xn.data(), MySim.pn.data(), MySim.xn.size());


			ImPlot::EndPlot();
		}
		ImGui::SeparatorText("Temperature Profile");
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);

		if (ImPlot::BeginPlot("Temperature Plot", ImVec2(-1, 500))) {


			ImPlot::PlotLine("T1 vs x Plot :", MySim.x1.data(), MySim.T1.data(), MySim.x1.size());
			ImPlot::PlotLine("T2 vs x Plot :", MySim.x2.data(), MySim.T2.data(), MySim.x2.size());
			ImPlot::PlotLine("T3 vs x Plot :", MySim.x3.data(), MySim.T3.data(), MySim.x3.size());
			ImPlot::PlotLine("T4 vs x Plot :", MySim.x4.data(), MySim.T4.data(), MySim.x4.size());
			ImPlot::PlotLine("T5 vs x Plot :", MySim.xn.data(), MySim.Tn.data(), MySim.xn.size());


			ImPlot::EndPlot();
		}
		ImGui::SeparatorText("Pressure compair");
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);

		if (ImPlot::BeginPlot("pressure compair Plot", ImVec2(-1, 500))) {


			ImPlot::PlotLine("p1/p2 plot,from fomular :", MySim.longx.data(), MySim.pp1.data(), MySim.longx.size());
			ImPlot::PlotLine("p1/p2 plot,from simulation :", MySim.longx.data(), MySim.ppp1.data(), MySim.longx.size());
			
			


			ImPlot::EndPlot();
		}
	
		

	}
	else
	{
		ImGui::SeparatorText("New Velocity Profile");
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);
		if (ImPlot::BeginPlot("new velocity Plot",ImVec2(-1,700))) 
		{


			ImPlot::PlotLine("velocity in magnetic field:", MySim.Veloci.data(), MySim.Z_axis.data(), MySim.Veloci.size());
			ImPlot::PlotLine("velocity in no magnetic field :", MySim.magnetic.data(), MySim.Z_axis.data(), MySim.magnetic.size());
			//ImPlot::PlotLine("magnetic field :", MySim.magnetic_fie.data(), MySim.Z_axis.data(), MySim.magnetic_fie.size());

			ImPlot::EndPlot();
		}
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);
		if (ImPlot::BeginPlot("magnetic Plot", ImVec2(-1, 500)))
		{

			ImPlot::PlotLine("magnetic field :", MySim.magnetic_fie.data(), MySim.Z_axis.data(), MySim.magnetic_fie.size());

			ImPlot::EndPlot();
		}
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);
		if (ImPlot::BeginPlot("Current Plot", ImVec2(-1, 500)))
		{

			ImPlot::PlotLine("Current :", MySim.current.data(), MySim.Z_axis.data(), MySim.current.size());

			ImPlot::EndPlot();
		}
		ImPlot::SetNextAxisToFit(ImAxis_X1);
		ImPlot::SetNextAxisToFit(ImAxis_Y1);
		if (ImPlot::BeginPlot("Electric Plot", ImVec2(-1, 500)))
		{

			ImPlot::PlotLine("Electric Field  :", MySim.Electric.data(), MySim.Z_axis.data(), MySim.Electric.size());

			ImPlot::EndPlot();
		}
	}
	
	ImGui::EndChild();
	ImGui::End();
}


