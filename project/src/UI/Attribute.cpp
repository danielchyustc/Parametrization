#include <UI/Attribute.h>

#include <UI/Grid.h>

#include <Engine/MeshEdit/Glue.h>
#include <Engine/MeshEdit/MinSurf.h>
#include <Engine/MeshEdit/IsotropicRemeshing.h>
#include <Engine/MeshEdit/Preprocessing.h>

#include <Engine/MeshEdit/Parameterization/Parameterization.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/PlanarEmbedding.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/TutteEmbedding.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/ASAP.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/ARAP.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/Progressive.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/Kinematic.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/VLG.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/CHM.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SphericalEmbedding.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SBE.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SCE.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVA.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVB.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVC.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SARAP.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVARAP.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/GaussMap.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/ProjTutte.h>
#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/ProjASAP.h>

#include <Engine/Scene/SObj.h>
#include <Engine/Scene/AllComponents.h>

#include <Engine/Primitive/Sphere.h>
#include <Engine/Primitive/Plane.h>
#include <Engine/Primitive/TriMesh.h>
#include <Engine/Primitive/Disk.h>
#include <Engine/Primitive/Capsule.h>

#include <Engine/Material/AllBSDFs.h>
#include <Engine/Material/Gooch.h>

#include <Engine/Light/AreaLight.h>
#include <Engine/Light/PointLight.h>
#include <Engine/Light/DirectionalLight.h>
#include <Engine/Light/SpotLight.h>
#include <Engine/Light/InfiniteAreaLight.h>
#include <Engine/Light/SphereLight.h>
#include <Engine/Light/DiskLight.h>
#include <Engine/Light/CapsuleLight.h>

#include <OpenGL/Camera.h>

#include <Basic/Math.h>
#include <Basic/EventManager.h>

#include <UDP/Visitor/MultiVisitor.h>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QLabel>
#include <Qt/RawAPI_OGLW.h>
#include <qcombobox.h>

#include <tuple>

using namespace Ubpa;

using namespace Ubpa::Math;

using namespace Ubpa;
using namespace std;

class Attribute::ComponentVisitor : public HeapObj,
	public SharedPtrMultiVisitor<Attribute::ComponentVisitor,
	Component, Light, Primitive, Material>
{
public:
	ComponentVisitor(Attribute* attr) : attr(attr) {
		Regist<CmptCamera,
			CmptGeometry,
			CmptLight,
			CmptMaterial,
			CmptTransform>();

		Regist<Sphere,
			Plane,
			TriMesh,
			Disk,
			Capsule>();

		Regist<AreaLight,
			PointLight,
			DirectionalLight,
			SpotLight,
			InfiniteAreaLight,
			SphereLight,
			DiskLight,
			CapsuleLight>();

		Regist<BSDF_Diffuse,
			BSDF_Glass,
			BSDF_Mirror,
			BSDF_Emission,
			BSDF_CookTorrance,
			BSDF_MetalWorkflow,
			BSDF_FrostedGlass,
			BSDF_Frostbite,
			Gooch>();
	}

public:
	static const Ptr<ComponentVisitor> New(Attribute* attr) {
		return Ubpa::New<ComponentVisitor>(attr);
	}

protected:
	virtual ~ComponentVisitor() = default;

protected:
	void ImplVisit(Ptr<CmptCamera> cmpt);

	void ImplVisit(Ptr<CmptGeometry> cmpt);
	void ImplVisit(Ptr<Sphere> sphere);
	void ImplVisit(Ptr<Plane> plane);
	void ImplVisit(Ptr<TriMesh> mesh);
	void ImplVisit(Ptr<Disk> disk);
	void ImplVisit(Ptr<Capsule> capsule);

	void ImplVisit(Ptr<CmptLight> cmpt);
	void ImplVisit(Ptr<AreaLight> light);
	void ImplVisit(Ptr<PointLight> light);
	void ImplVisit(Ptr<DirectionalLight> light);
	void ImplVisit(Ptr<SpotLight> light);
	void ImplVisit(Ptr<InfiniteAreaLight> light);
	void ImplVisit(Ptr<SphereLight> light);
	void ImplVisit(Ptr<DiskLight> light);
	void ImplVisit(Ptr<CapsuleLight> light);

	void ImplVisit(Ptr<CmptMaterial> cmpt);
	void ImplVisit(Ptr<BSDF_Diffuse> bsdf);
	void ImplVisit(Ptr<BSDF_Emission> bsdf);
	void ImplVisit(Ptr<BSDF_Glass> bsdf);
	void ImplVisit(Ptr<BSDF_Mirror> bsdf);
	void ImplVisit(Ptr<BSDF_CookTorrance> bsdf);
	void ImplVisit(Ptr<BSDF_MetalWorkflow> bsdf);
	void ImplVisit(Ptr<BSDF_FrostedGlass> bsdf);
	void ImplVisit(Ptr<Gooch> gooch);
	void ImplVisit(Ptr<BSDF_Frostbite> bsdf);

	void ImplVisit(Ptr<CmptTransform> cmpt);

private:
	QWidget* GenItem(Ptr<Component> component, const QString& str) {
		auto item = new QWidget;
		attr->componentType2item[typeid(*component)] = item;
		attr->tbox->insertItem(0, item, str);
		return item;
	}

	Ptr<Grid> GetGrid(QWidget* item) {
		auto target = attr->item2grid.find(item);
		if (target != attr->item2grid.end())
			return target->second;

		auto grid = Grid::New(item);
		attr->item2grid[item] = grid;

		return grid;
	}

private:
	Attribute* attr;
};

// -------------- Transform --------------

void Attribute::ComponentVisitor::ImplVisit(Ptr<CmptTransform> transform) {
	auto item = GenItem(transform, "Transform");
	auto grid = GetGrid(item);

	// position
	grid->AddText("- Position");
	grid->AddEditVal({ "x","y","z" }, transform->GetPosition().cast_to<Ubpa::valf3>(), Ubpa::valf3(0.1f), [=](const Ubpa::valf3& val) {
		transform->SetPosition(val.cast_to<Ubpa::pointf3>());
		});
	grid->AddText("- Rotation");
	auto e = transform->GetRotationEuler().cast_to<Ubpa::vecf3>() / Ubpa::PI<float> * 180;
	grid->AddEditVal({ "x","y","z" }, e.cast_to<Ubpa::valf3>(), Ubpa::valf3(-90, -180, -180), Ubpa::valf3(90, 180, 180), Ubpa::vali3(180, 360, 360), [=](const Ubpa::valf3& val) {
		transform->SetRotation(Ubpa::eulerf{ Ubpa::to_radian(val[0]),Ubpa::to_radian(val[1]) ,Ubpa::to_radian(val[2]) });
		});

	grid->AddText("- Scale");
	grid->AddEditVal({ "x","y","z" }, transform->GetScale().cast_to<Ubpa::valf3>(), Ubpa::valf3(0.1f), [=](const Ubpa::valf3& val) {
		transform->SetScale(val.cast_to<Ubpa::scalef3>());
		});
}

// -------------- Camera --------------

void Attribute::ComponentVisitor::ImplVisit(Ptr<CmptCamera> camera) {
	auto item = GenItem(camera, "Camera");
	auto grid = GetGrid(item);

	grid->AddEditVal("- Field of View", camera->GetFOV(), 1, 179, [camera](double fov) {
		camera->SetFOV(fov);
		});

	grid->AddButton("- To Roamer Camera", [camera]() {
		WPtr<Camera> wptrRoamerCam;
		if (!GS::GetV("Roamer::camera", wptrRoamerCam)) {
			printf("ERROR::Attribute::ComponentVisitor::ImplVisit(Ptr<CmptCamera>):\n"
				"\t""GS::GetV(\"Roamer::camera\") fail\n");
			return;
		}

		auto roamerCam = wptrRoamerCam.lock();
		if (!roamerCam) {
			printf("ERROR::Attribute::ComponentVisitor::ImplVisit(Ptr<CmptCamera>):\n"
				"\t""roamerCam is nullptr\n");
			return;
		}

		auto fov = camera->GetFOV();
		auto transform = camera->GetSObj()->GetLocalToWorldMatrix();

		auto eulerAngle = transform.decompose_euler();
		roamerCam->SetPose(transform.decompose_position(), Ubpa::to_radian(eulerAngle[1]), Ubpa::to_radian(eulerAngle[0]));
		roamerCam->SetFOV(fov);
		EventMngr::GetInstance().Response(0, roamerCam);
		});

	grid->AddButton("- From Roamer Camera", [camera]() {
		WPtr<Camera> wptrRoamerCam;
		if (!GS::GetV("Roamer::camera", wptrRoamerCam)) {
			printf("ERROR::Attribute::ComponentVisitor::ImplVisit(Ptr<CmptCamera>):\n"
				"\t""GS::GetV(\"Roamer::camera\") fail\n");
			return;
		}

		auto roamerCam = wptrRoamerCam.lock();
		if (!roamerCam) {
			printf("ERROR::Attribute::ComponentVisitor::ImplVisit(Ptr<CmptCamera>):\n"
				"\t""roamerCam is nullptr\n");
			return;
		}

		auto transformCmpt = camera->GetSObj()->GetComponent<CmptTransform>();
		if (!transformCmpt)
			transformCmpt = CmptTransform::New(camera->GetSObj());

		camera->SetFOV(roamerCam->GetFOV());

		auto front = roamerCam->GetFront();
		auto pos = roamerCam->GetPos();
		auto lookAt = Ubpa::transformf::look_at(pos, pos + front);
		auto w2parent = camera->GetSObj()->GetLocalToWorldMatrix().inverse() * transformCmpt->GetTransform();
		transformCmpt->SetTransform(lookAt.inverse() * w2parent);

		});
}

// -------------- Geometry --------------

void Attribute::ComponentVisitor::ImplVisit(Ptr<CmptGeometry> geo) {
	auto item = GenItem(geo, "Geometry");
	auto grid = GetGrid(item);

	string typeStr;
	BasicSharedPtrVisitor<Primitive> getTypeStr;
	getTypeStr.Regist(
		[&typeStr](Ptr<Sphere>) { typeStr = "Sphere"; },
		[&typeStr](Ptr<Plane>) { typeStr = "Plane"; },
		[&typeStr](Ptr<TriMesh>) { typeStr = "TriMesh"; },
		[&typeStr](Ptr<Disk>) { typeStr = "Disk"; },
		[&typeStr](Ptr<Capsule>) { typeStr = "Capsule"; }
	);

	Grid::pSlotMap pSlotMap(new Grid::SlotMap);
	Grid::wpSlotMap wpSlotMap = pSlotMap;

	QComboBox* combobox = new QComboBox;

	(*pSlotMap)["None"] = [=]() {
		grid->Clear();
		grid->AddComboBox("Type", "None", wpSlotMap.lock());
		geo->primitive = nullptr;
	};
	(*pSlotMap)["Sphere"] = [=]() {
		grid->Clear();
		grid->AddComboBox("Type", "Sphere", wpSlotMap.lock());

		auto sphere = Sphere::New();
		geo->primitive = sphere;
		Visit(geo->primitive);
	};
	(*pSlotMap)["Plane"] = [=]() {
		grid->Clear();
		grid->AddComboBox("Type", "Plane", wpSlotMap.lock());

		auto plane = Plane::New();
		geo->primitive = plane;
		Visit(geo->primitive);
	};
	(*pSlotMap)["TriMesh"] = [=]() {
		// not support now
		grid->Clear();
		grid->AddComboBox("Type", "None", wpSlotMap.lock());
		geo->primitive = nullptr;
	};
	(*pSlotMap)["Disk"] = [=]() {
		grid->Clear();
		grid->AddComboBox("Type", "Disk", wpSlotMap.lock());

		auto disk = Disk::New();
		geo->primitive = disk;
		Visit(geo->primitive);
	};
	(*pSlotMap)["Capsule"] = [=]() {
		grid->Clear();
		grid->AddComboBox("Type", "Capsule", wpSlotMap.lock());

		auto capsule = Capsule::New();
		geo->primitive = capsule;
		Visit(geo->primitive);
	};

	if (geo->primitive) {
		getTypeStr.Visit(geo->primitive);
		grid->AddComboBox(combobox, "Type", typeStr, pSlotMap);

		Visit(geo->primitive);
	}
	else
		grid->AddComboBox(combobox, "Type", "None", pSlotMap);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<Sphere> sphere) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptGeometry)]);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<Plane> plane) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptGeometry)]);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<Disk> disk) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptGeometry)]);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<TriMesh> mesh) {
	auto preproc = Preprocessing::New(mesh);

	auto grid = GetGrid(attr->componentType2item[typeid(CmptGeometry)]);
	grid->AddText("- Triangle", mesh->GetIndice().size() / 3);
	grid->AddText("- Vertex", mesh->GetPositions().size());

	if (preproc->IsClosed())
	{
		grid->AddTitle("- Spherical Embedding");

		Grid::pSlotMap sphere_slotmap = std::make_shared<Grid::SlotMap>();
		(*sphere_slotmap)["Gauss Map"] = [this]() { attr->sphere_method_ = kGaussMap; };
		(*sphere_slotmap)["ProjTutte"] = [this]() { attr->sphere_method_ = kProjTutte; };
		(*sphere_slotmap)["ProjASAP"] = [this]() { attr->sphere_method_ = kProjASAP; };
		(*sphere_slotmap)["SBE"] = [this]() { attr->sphere_method_ = kSBE; };
		(*sphere_slotmap)["SCE"] = [this]() { attr->sphere_method_ = kSCE; };
		(*sphere_slotmap)["CdVA"] = [this]() { attr->sphere_method_ = kCdVA; };
		(*sphere_slotmap)["CdVB"] = [this]() { attr->sphere_method_ = kCdVB; };
		(*sphere_slotmap)["CdVC"] = [this]() { attr->sphere_method_ = kCdVC; };
		(*sphere_slotmap)["CdVARAP"] = [this]() { attr->sphere_method_ = kCdVARAP; };
		(*sphere_slotmap)["SARAP"] = [this]() { attr->sphere_method_ = kSARAP; };
		grid->AddComboBox("Method", "Gauss Map", sphere_slotmap);
		attr->sphere_method_ = kGaussMap;

		Grid::pSlotMap sphere_view_mode_slotmap = std::make_shared<Grid::SlotMap>();
		(*sphere_view_mode_slotmap)["Mesh"] = [=]() {
			attr->view_mode_ = kMesh;
			if (attr->param_ != nullptr)
				attr->param_->SetViewMode(kMesh);
			attr->pOGLW->DirtyVAO(mesh);
		};
		(*sphere_view_mode_slotmap)["Sphere"] = [=]() {
			attr->view_mode_ = kSphere;
			if (attr->param_ != nullptr)
				attr->param_->SetViewMode(kSphere);
			attr->pOGLW->DirtyVAO(mesh);
		};
		(*sphere_view_mode_slotmap)["Plane"] = [=]() {
			attr->view_mode_ = kPlane;
			if (attr->param_ != nullptr)
				attr->param_->SetViewMode(kPlane);
			attr->pOGLW->DirtyVAO(mesh);
		};
		grid->AddComboBox("View Mode", "Mesh", sphere_view_mode_slotmap);

		grid->AddButton("Initialize", [=]() {
			switch (attr->sphere_method_)
			{
			case kGaussMap:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = GaussMap::New(preproc);
				else
					attr->param_ = GaussMap::New(attr->param_);
				break;
			case kProjTutte:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = ProjTutte::New(preproc);
				else
					attr->param_ = ProjTutte::New(attr->param_);
				break;
			case kProjASAP:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = ProjASAP::New(preproc);
				else
					attr->param_ = ProjASAP::New(attr->param_);
				break;
			case kSBE:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = SBE::New(preproc);
				else
					attr->param_ = SBE::New(attr->param_);
				break;
			case kSCE:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = SCE::New(preproc);
				else
					attr->param_ = SCE::New(attr->param_);
				break;
			case kCdVA:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = CdVA::New(preproc);
				else
					attr->param_ = CdVA::New(attr->param_);
				break;
			case kCdVB:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = CdVB::New(preproc);
				else
					attr->param_ = CdVB::New(attr->param_);
				break;
			case kCdVC:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = CdVC::New(preproc);
				else
					attr->param_ = CdVC::New(attr->param_);
				break;
			case kCdVARAP:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = CdVARAP::New(preproc);
				else
					attr->param_ = CdVARAP::New(attr->param_);
				break;
			case kSARAP:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = SARAP::New(preproc);
				else
					attr->param_ = SARAP::New(attr->param_);
				break;
			default:
				break;
			}
			printf("Gauss Map\n");
			attr->param_->SetViewMode(attr->view_mode_);
			printf("View mode set.\n");
			if (attr->param_->Run())
				printf("Spherical embedding initialization done\n");
			attr->pOGLW->DirtyVAO(mesh);
			});

		grid->AddButton("Iterate", [=]() {
			if (attr->param_ == nullptr)
			{
				printf("hasn't initialize spherical embedding yet!\n");
				return;
			}
			else if (attr->param_->Iterate())
				printf("Spheriical Embedding Iteration done\n");
			attr->pOGLW->DirtyVAO(mesh);
			});

		grid->AddButton("Iterate Till Done", [=]() {
			if (attr->param_ == nullptr)
			{
				printf("hasn't initialize spherical embedding yet!\n");
				return;
			}
			else if (attr->param_->IterateTillDone())
				printf("Spherical Embedding Iteration finished\n");
			attr->pOGLW->DirtyVAO(mesh);
			});

		grid->AddButton("Move Barycenter to Origin", [=]() {
			if (attr->param_ == nullptr)
			{
				printf("hasn't initialize spherical embedding yet!\n");
				return;
			}
			attr->param_->MoveBarycenterToOrigin();
			attr->param_->UpdateTriMesh();
			attr->pOGLW->DirtyVAO(mesh);
			});

		grid->AddButton("Area Correction", [=]() {
			if (attr->param_ == nullptr)
			{
				printf("hasn't initialize spherical embedding yet!\n");
				return;
			}
			attr->param_->AreaCorrection();
			attr->param_->UpdateTriMesh();
			attr->pOGLW->DirtyVAO(mesh);
			});

		vector<string> weight_items{ "Equal", "Cotangent" };
		vector<string> shape_items{ "Circle", "Square" };


		auto setWeightSlot = [=](const string& item) {
			if (item == weight_items[0])
			{
				attr->SetWeightType(kEqual);
			}
			else if (item == weight_items[1])
			{
				attr->SetWeightType(kCotangent);
			}
		};

		auto setShapeSlot = [=](const string& item) {
			if (item == "Circle")
			{
				attr->SetBoundShape(kCircle);
			}
			else if (item == "Square")
			{
				attr->SetBoundShape(kSquare);
			}
		};

		grid->AddComboBox("- Weight Type", weight_items, "Change Weight", setWeightSlot);
	}
	else 
	{
		grid->AddButton("Minimize Surface", [=]() {
			auto minSurf = MinSurf::New(mesh, attr->weight_type_);
			minSurf->Run();
			attr->pOGLW->DirtyVAO(mesh);
			});

		grid->AddButton("Isotropic Remeshing", [mesh, pOGLW = attr->pOGLW]() {
			printf("[Isotropic Remeshing] start\n");
			auto isotropicRemeshing = IsotropicRemeshing::New(mesh);
			printf("[Isotropic Remeshing] init done\n");
			if (isotropicRemeshing->Run(5))
				printf("[Isotropic Remeshing] success\n");
			else
				printf("[Isotropic Remeshing] fail\n");
			pOGLW->DirtyVAO(mesh);
		});

		grid->AddTitle("[ Parametrization ]");
		grid->AddTitle("- Planar Embedding");

		Grid::pSlotMap plane_slotmap = std::make_shared<Grid::SlotMap>();
		(*plane_slotmap)["Tutte Embedding"] = [this]() { attr->plane_method_ = kTutteEmbedding; };
		(*plane_slotmap)["ASAP"] = [this]() { attr->plane_method_ = kASAP; };
		(*plane_slotmap)["ARAP"] = [this]() { attr->plane_method_ = kARAP; };
		(*plane_slotmap)["Progressive"] = [this]() { attr->plane_method_ = kProg; };
		(*plane_slotmap)["Kinematic"] = [this]() { attr->plane_method_ = kKine; };
		(*plane_slotmap)["VLG"] = [this]() { attr->plane_method_ = kVLG; };
		(*plane_slotmap)["CHM"] = [this]() { attr->plane_method_ = kCHM; };
		grid->AddComboBox("Method", "Tutte Embedding", plane_slotmap);
		attr->plane_method_ = kTutteEmbedding;

		Grid::pSlotMap plane_view_mode_slotmap = std::make_shared<Grid::SlotMap>();
		(*plane_view_mode_slotmap)["Mesh"] = [=]() {
			attr->view_mode_ = kMesh;
			if (attr->param_ != nullptr)
				attr->param_->SetViewMode(kMesh);
			attr->pOGLW->DirtyVAO(mesh);
		};
		(*plane_view_mode_slotmap)["Plane"] = [=]() {
			attr->view_mode_ = kPlane;
			if (attr->param_ != nullptr)
				attr->param_->SetViewMode(kPlane);
			attr->pOGLW->DirtyVAO(mesh);
		};
		grid->AddComboBox("View Mode", "Mesh", plane_view_mode_slotmap);

		grid->AddButton("Initialize", [=]() {
			switch (attr->plane_method_)
			{
			case kTutteEmbedding:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = TutteEmbedding::New(preproc);
				else
					attr->param_ = TutteEmbedding::New(attr->param_);
				break;
			case kASAP:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = ASAP::New(preproc);
				else
					attr->param_ = ASAP::New(attr->param_);
				break;
			case kARAP:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
				{
					attr->param_ = ASAP::New(preproc);
					attr->param_->SetViewMode(attr->view_mode_);
					if (attr->param_->Run())
						printf("Planar embedding initialization done\n");
					attr->pOGLW->DirtyVAO(mesh);
					attr->param_->Run();
				}
				attr->param_ = ARAP::New(attr->param_);
				grid->AddButton("Iterate", [=]() {
					if (attr->param_ == nullptr)
					{
						printf("hasn't initialize planar embedding yet!\n");
						return;
					}
					else if (attr->param_->Iterate())
						printf("Planar Embedding Iteration done\n");
					attr->pOGLW->DirtyVAO(mesh);
					});
				break;
			case kProg:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = Progressive::New(preproc);
				else
					attr->param_ = Progressive::New(attr->param_);
				grid->AddButton("Iterate", [=]() {
					if (attr->param_ == nullptr)
					{
						printf("hasn't initialize planar embedding yet!\n");
						return;
					}
					else if (attr->param_->Iterate())
						printf("Planar Embedding Iteration done\n");
					attr->pOGLW->DirtyVAO(mesh);
					});
				break;
			case kKine:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = Kinematic::New(preproc);
				else
					attr->param_ = Kinematic::New(attr->param_);
				grid->AddButton("Iterate", [=]() {
					if (attr->param_ == nullptr)
					{
						printf("hasn't initialize planar embedding yet!\n");
						return;
					}
					else if (attr->param_->Iterate())
						printf("Planar Embedding Iteration done\n");
					attr->pOGLW->DirtyVAO(mesh);
					});
				break;
			case kVLG:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = VLG::New(preproc);
				else
					attr->param_ = VLG::New(attr->param_);
				grid->AddButton("Iterate", [=]() {
					if (attr->param_ == nullptr)
					{
						printf("hasn't initialize planar embedding yet!\n");
						return;
					}
					else if (attr->param_->Iterate())
						printf("Planar Embedding Iteration done\n");
					attr->pOGLW->DirtyVAO(mesh);
					});
				break;
			case kCHM:
				if (attr->param_ == nullptr || attr->param_->GetTriMesh() != mesh)
					attr->param_ = CHM::New(preproc);
				else
					attr->param_ = CHM::New(attr->param_);
				grid->AddButton("Iterate", [=]() {
					if (attr->param_ == nullptr)
					{
						printf("hasn't initialize planar embedding yet!\n");
						return;
					}
					else if (attr->param_->Iterate())
						printf("Planar Embedding Iteration done\n");
					attr->pOGLW->DirtyVAO(mesh);
					});
				break;
			default:
				break;
			}
			attr->param_->SetViewMode(attr->view_mode_);
			if (attr->param_->Run())
				printf("Planar embedding initialization done\n");
			attr->pOGLW->DirtyVAO(mesh);
			});

		vector<string> weight_items{ "Equal", "Cotangent" };
		vector<string> shape_items{ "Circle", "Square" };


		auto setWeightSlot = [=](const string& item) {
			if (item == weight_items[0])
			{
				attr->SetWeightType(kEqual);
			}
			else if (item == weight_items[1])
			{
				attr->SetWeightType(kCotangent);
			}
		};

		auto setShapeSlot = [=](const string& item) {
			if (item == "Circle")
			{
				attr->SetBoundShape(kCircle);
			}
			else if (item == "Square")
			{
				attr->SetBoundShape(kSquare);
			}
		};

		grid->AddComboBox("- Weight Type", weight_items, "Change Weight", setWeightSlot);
		grid->AddComboBox("- Bound Shape", shape_items, "Change Shape", setShapeSlot);
	}

}

void Attribute::ComponentVisitor::ImplVisit(Ptr<Capsule> capsule) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptGeometry)]);
	grid->AddEditVal("- Height", capsule->height, 0.f, 100.f, 1000);
}

// -------------- Material --------------

void Attribute::ComponentVisitor::ImplVisit(Ptr<CmptMaterial> cmpt) {
	auto item = GenItem(cmpt, "Material");
	auto grid = GetGrid(item);

	string typeStr;
	BasicSharedPtrVisitor<Material> getTypeStr;
	getTypeStr.Regist(
		[&typeStr](Ptr<BSDF_Diffuse>) { typeStr = "BSDF_Diffuse"; },
		[&typeStr](Ptr<BSDF_Emission>) { typeStr = "BSDF_Emission"; },
		[&typeStr](Ptr<BSDF_Glass>) { typeStr = "BSDF_Glass"; },
		[&typeStr](Ptr<BSDF_Mirror>) { typeStr = "BSDF_Mirror"; },
		[&typeStr](Ptr<BSDF_CookTorrance>) { typeStr = "BSDF_CookTorrance"; },
		[&typeStr](Ptr<BSDF_MetalWorkflow>) { typeStr = "BSDF_MetalWorkflow"; },
		[&typeStr](Ptr<BSDF_FrostedGlass>) { typeStr = "BSDF_FrostedGlass"; },
		[&typeStr](Ptr<Gooch>) { typeStr = "Gooch"; },
		[&typeStr](Ptr<BSDF_Frostbite>) { typeStr = "BSDF_Frostbite"; }
	);

	const int materialNum = 10;
	tuple<string, function<Ptr<Material>()>> bsdfArr[materialNum] = {
		{"None", []()->Ptr<Material> { return nullptr; } },
		{"BSDF_Diffuse", []()->Ptr<Material> { return BSDF_Diffuse::New(); } },
		{"BSDF_Emission", []()->Ptr<Material> { return BSDF_Emission::New(); } },
		{"BSDF_Glass", []()->Ptr<Material> { return BSDF_Glass::New(); } },
		{"BSDF_Mirror", []()->Ptr<Material> { return BSDF_Mirror::New(); } },
		{"BSDF_CookTorrance", []()->Ptr<Material> { return BSDF_CookTorrance::New(); } },
		{"BSDF_MetalWorkflow", []()->Ptr<Material> { return BSDF_MetalWorkflow::New(); } },
		{"BSDF_FrostedGlass", []()->Ptr<Material> { return BSDF_FrostedGlass::New(); } },
		{"Gooch", []()->Ptr<Material> { return Gooch::New(); } },
		{"BSDF_Frostbite", []()->Ptr<Material> { return BSDF_Frostbite::New(); } },
	};

	Grid::pSlotMap pSlotMap(new Grid::SlotMap);
	Grid::wpSlotMap wpSlotMap = pSlotMap;

	for (int i = 0; i < materialNum; i++) {
		(*pSlotMap)[get<0>(bsdfArr[i])] = [=]() {
			grid->Clear();
			grid->AddComboBox("Type", get<0>(bsdfArr[i]), wpSlotMap.lock());

			auto bsdf = get<1>(bsdfArr[i])();
			cmpt->material = bsdf;
			if (bsdf != nullptr)
				Visit(bsdf);
		};
	}

	if (cmpt->material) {
		getTypeStr.Visit(cmpt->material);
		grid->AddComboBox("Type", typeStr, pSlotMap);

		Visit(cmpt->material);
	}
	else
		grid->AddComboBox("Type", "None", pSlotMap);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_Diffuse> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);
	grid->AddEditColor("- Albedo", bsdf->colorFactor);
	grid->AddEditImage("- Albedo Texture", bsdf->albedoTexture);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_Emission> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);
	grid->AddEditColor("- Color", bsdf->color);
	grid->AddEditVal("- Intensity", bsdf->intensity, 0, 100, 1000);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_Glass> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);
	grid->AddEditColor("- Transmittance", bsdf->transmittance);
	grid->AddEditColor("- Reflectance", bsdf->reflectance);
	grid->AddEditVal("- ior", bsdf->ior, 0.01);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_Mirror> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);
	grid->AddEditColor("- Reflectance", bsdf->reflectance);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_CookTorrance> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);
	grid->AddEditColor("- Reflectance", bsdf->refletance);
	grid->AddEditColor("- Albedo", bsdf->albedo);
	grid->AddEditVal("- Index of Refract", bsdf->ior, 0.01);
	grid->AddEditVal("- Roughness", bsdf->m, 0, 1, 100);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_MetalWorkflow> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);

	grid->AddEditColor("- Albedo Color", bsdf->colorFactor);
	grid->AddEditImage("- Albedo Texture", bsdf->albedoTexture);

	grid->AddEditVal("- Metallic Factor", bsdf->metallicFactor, 0, 1, 100);
	grid->AddEditImage("- Metallic Texture", bsdf->metallicTexture);

	grid->AddEditVal("- Roughness Factor", bsdf->roughnessFactor, 0, 1, 100);
	grid->AddEditImage("- Roughness Texture", bsdf->roughnessTexture);

	grid->AddEditImage("- AO Texture", bsdf->aoTexture);

	grid->AddEditImage("- Normal Texture", bsdf->normalTexture);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_FrostedGlass> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);

	grid->AddEditVal("- IOR", bsdf->ior, 1., 20., 190);

	grid->AddEditColor("- Color Factor", bsdf->colorFactor);
	grid->AddEditImage("- Color Texture", bsdf->colorTexture);

	grid->AddEditVal("- Roughness Factor", bsdf->roughnessFactor, 0.02, 1., 98);
	grid->AddEditImage("- Roughness Texture", bsdf->roughnessTexture);

	grid->AddEditImage("- AO Texture", bsdf->aoTexture);
	grid->AddEditImage("- Normal Texture", bsdf->normalTexture);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<Gooch> gooch) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);

	grid->AddEditColor("- Color Factor", gooch->colorFactor);
	grid->AddEditImage("- Color Texture", gooch->colorTexture);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<BSDF_Frostbite> bsdf) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptMaterial)]);

	grid->AddEditColor("- Albedo Color", bsdf->colorFactor);
	grid->AddEditImage("- Albedo Texture", bsdf->albedoTexture);

	grid->AddEditVal("- Metallic Factor", bsdf->metallicFactor, 0, 1, 100);
	grid->AddEditImage("- Metallic Texture", bsdf->metallicTexture);

	grid->AddEditVal("- Roughness Factor", bsdf->roughnessFactor, 0, 1, 100);
	grid->AddEditImage("- Roughness Texture", bsdf->roughnessTexture);

	grid->AddEditImage("- AO Texture", bsdf->aoTexture);

	grid->AddEditImage("- Normal Texture", bsdf->normalTexture);
}

// -------------- Light --------------

void Attribute::ComponentVisitor::ImplVisit(Ptr<CmptLight> cmpt) {
	auto item = GenItem(cmpt, "Light");
	auto grid = GetGrid(item);

	string typeStr;
	BasicSharedPtrVisitor<Light> getTypeStr;
	getTypeStr.Regist(
		[&typeStr](Ptr<AreaLight>) { typeStr = "AreaLight"; },
		[&typeStr](Ptr<PointLight>) { typeStr = "PointLight"; },
		[&typeStr](Ptr<DirectionalLight>) { typeStr = "DirectionalLight"; },
		[&typeStr](Ptr<SpotLight>) { typeStr = "SpotLight"; },
		[&typeStr](Ptr<InfiniteAreaLight>) { typeStr = "InfiniteAreaLight"; },
		[&typeStr](Ptr<SphereLight>) { typeStr = "SphereLight"; },
		[&typeStr](Ptr<DiskLight>) { typeStr = "DiskLight"; },
		[&typeStr](Ptr<CapsuleLight>) { typeStr = "CapsuleLight"; }
	);

	const int lightNum = 9;
	tuple<string, function<Ptr<Light>()>> lightArr[lightNum] = {
		{"None", []()->Ptr<Light> { return nullptr; } },
		{"AreaLight", []()->Ptr<Light> { return AreaLight::New(); } },
		{"PointLight", []()->Ptr<Light> { return PointLight::New(); } },
		{"DirectionalLight", []()->Ptr<Light> { return DirectionalLight::New(); } },
		{"SpotLight", []()->Ptr<Light> { return SpotLight::New(); } },
		{"InfiniteAreaLight", []()->Ptr<Light> { return InfiniteAreaLight::New(nullptr); } },
		{"SphereLight", []()->Ptr<Light> { return SphereLight::New(); } },
		{"DiskLight", []()->Ptr<Light> { return DiskLight::New(); } },
		{"CapsuleLight", []()->Ptr<Light> { return CapsuleLight::New(); } },
	};

	Grid::pSlotMap pSlotMap(new Grid::SlotMap);
	Grid::wpSlotMap wpSlotMap = pSlotMap;

	for (int i = 0; i < lightNum; i++) {
		(*pSlotMap)[get<0>(lightArr[i])] = [=]() {
			grid->Clear();
			grid->AddComboBox("Type", get<0>(lightArr[i]), wpSlotMap.lock());

			auto lightbase = get<1>(lightArr[i])();
			cmpt->light = lightbase;
			if (lightbase != nullptr)
				Visit(lightbase);
		};
	}

	if (cmpt->light) {
		getTypeStr.Visit(cmpt->light);
		grid->AddComboBox("Type", typeStr, pSlotMap);

		Visit(cmpt->light);
	}
	else
		grid->AddComboBox("Type", "None", pSlotMap);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<AreaLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
	grid->AddEditVal("- Width", light->width, 0.1);
	grid->AddEditVal("- Height", light->height, 0.1);
	grid->AddButton("- Autocorrection", [this, light]() {
		auto curSObj = attr->GetCurSObj();
		if (!curSObj)
			return;

		// geo
		auto geoCmpt = curSObj->GetComponent<CmptGeometry>();
		if (!geoCmpt)
			geoCmpt = CmptGeometry::New(curSObj, Plane::New());
		if (!CastTo<Plane>(geoCmpt->primitive))
			geoCmpt->primitive = Plane::New();

		// transform
		auto tsfmCmpt = curSObj->GetComponent<CmptTransform>();
		if (!tsfmCmpt)
			tsfmCmpt = CmptTransform::New(curSObj);
		tsfmCmpt->SetScale({ light->width, 1.f, light->height });

		// material
		auto matCmpt = curSObj->GetComponent<CmptMaterial>();
		if (!matCmpt)
			matCmpt = CmptMaterial::New(curSObj, BSDF_Emission::New());
		auto emission = CastTo<BSDF_Emission>(matCmpt->material);
		if (!emission)
			matCmpt->material = emission = BSDF_Emission::New();
		emission->color = light->color;
		emission->intensity = light->Luminance().illumination() / light->color.illumination();

		attr->SetSObj(curSObj);
		attr->SetCurCmpt<CmptLight>();
		});
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<PointLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
	grid->AddEditVal("- Radius", light->radius, 0, 100.0, 1000);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<DirectionalLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<SpotLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
	grid->AddEditVal("- Radius", light->radius, 0, 100.0, 1000);
	grid->AddEditVal("- Angle", light->angle, 1.0, 179.0, 178);
	grid->AddEditVal("- Full Ratio", light->fullRatio, 0.0, 1.0, 100);
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<InfiniteAreaLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditVal("- Intensity", light->intensity, 0, 20, 2000);
	grid->AddEditColor("- Color Factor", light->colorFactor);
	grid->AddEditImage("- Envirment Image", light->GetImg(), [light](Ptr<Image> img) {
		light->SetImg(img);
		});
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<SphereLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Radius", light->radius, 0.01, 100.0, 1000);

	grid->AddButton("- Autocorrection", [this, light]() {
		auto curSObj = attr->GetCurSObj();
		if (!curSObj)
			return;

		// geo
		auto geoCmpt = curSObj->GetComponent<CmptGeometry>();
		if (!geoCmpt)
			geoCmpt = CmptGeometry::New(curSObj, Sphere::New());
		if (!CastTo<Sphere>(geoCmpt->primitive))
			geoCmpt->primitive = Sphere::New();

		// transform
		auto tsfmCmpt = curSObj->GetComponent<CmptTransform>();
		if (!tsfmCmpt)
			tsfmCmpt = CmptTransform::New(curSObj);
		tsfmCmpt->SetScale(light->radius);

		// material
		auto matCmpt = curSObj->GetComponent<CmptMaterial>();
		if (!matCmpt)
			matCmpt = CmptMaterial::New(curSObj, BSDF_Emission::New());
		auto emission = CastTo<BSDF_Emission>(matCmpt->material);
		if (!emission)
			matCmpt->material = emission = BSDF_Emission::New();
		emission->color = light->color;
		emission->intensity = light->Luminance().illumination() / light->color.illumination();

		attr->SetSObj(curSObj);
		attr->SetCurCmpt<CmptLight>();
		});
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<DiskLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Radius", light->radius, 0, 100, 1000);

	grid->AddButton("- Autocorrection", [this, light]() {
		auto curSObj = attr->GetCurSObj();
		if (!curSObj)
			return;

		// geo
		auto geoCmpt = curSObj->GetComponent<CmptGeometry>();
		if (!geoCmpt)
			geoCmpt = CmptGeometry::New(curSObj, Disk::New());
		if (!CastTo<Disk>(geoCmpt->primitive))
			geoCmpt->primitive = Disk::New();

		// transform
		auto tsfmCmpt = curSObj->GetComponent<CmptTransform>();
		if (!tsfmCmpt)
			tsfmCmpt = CmptTransform::New(curSObj);
		tsfmCmpt->SetScale({ light->radius, 1.0f, light->radius });

		// material
		auto matCmpt = curSObj->GetComponent<CmptMaterial>();
		if (!matCmpt)
			matCmpt = CmptMaterial::New(curSObj, BSDF_Emission::New());
		auto emission = CastTo<BSDF_Emission>(matCmpt->material);
		if (!emission)
			matCmpt->material = emission = BSDF_Emission::New();
		emission->color = light->color;
		emission->intensity = light->Luminance().illumination() / light->color.illumination();

		attr->SetSObj(curSObj);
		attr->SetCurCmpt<CmptLight>();
		});
}

void Attribute::ComponentVisitor::ImplVisit(Ptr<CapsuleLight> light) {
	auto grid = GetGrid(attr->componentType2item[typeid(CmptLight)]);
	grid->AddEditVal("- Intensity", light->intensity, 0, 10000, 100000);
	grid->AddEditColor("- Color", light->color);
	grid->AddEditVal("- Radius", light->radius, 0, 100, 1000);
	grid->AddEditVal("- Height", light->height, 0, 100, 1000);

	grid->AddButton("- Autocorrection", [this, light]() {
		auto curSObj = attr->GetCurSObj();
		if (!curSObj)
			return;

		// geo
		auto geoCmpt = curSObj->GetComponent<CmptGeometry>();
		if (!geoCmpt)
			geoCmpt = CmptGeometry::New(curSObj, Capsule::New());
		auto capsule = CastTo<Capsule>(geoCmpt->primitive);
		if (!capsule)
			geoCmpt->primitive = capsule = Capsule::New();
		capsule->height = light->height / light->radius;

		// transform
		auto tsfmCmpt = curSObj->GetComponent<CmptTransform>();
		if (!tsfmCmpt)
			tsfmCmpt = CmptTransform::New(curSObj);
		tsfmCmpt->SetScale(light->radius);

		// material
		auto matCmpt = curSObj->GetComponent<CmptMaterial>();
		if (!matCmpt)
			matCmpt = CmptMaterial::New(curSObj, BSDF_Emission::New());
		auto emission = CastTo<BSDF_Emission>(matCmpt->material);
		if (!emission)
			matCmpt->material = emission = BSDF_Emission::New();
		emission->color = light->color;
		emission->intensity = light->Luminance().illumination() / light->color.illumination();

		attr->SetSObj(curSObj);
		attr->SetCurCmpt<CmptLight>();
		});
}

// -------------- Attribute --------------

Attribute::Attribute()
	: tbox(nullptr), visitor(ComponentVisitor::New(this)) {
}

void Attribute::Init(QToolBox* tbox, RawAPI_OGLW* pOGLW) {
	this->tbox = tbox;
	this->pOGLW = pOGLW;
	SetSObj(nullptr);
}

void Attribute::SetSObj(Ptr<SObj> sobj) {
	if (!tbox)
		return;

	curSObj = sobj;

	// clear
	componentType2item.clear();

	for (auto pair : item2grid)
		pair.second->Clear();
	item2grid.clear();

	// 一个 component 对应一个 widget
	int num = tbox->count();
	while (num-- > 0) {
		// remove item 不会删除 item
		auto item = tbox->widget(0);
		tbox->removeItem(0);
		item->deleteLater();
	}

	if (sobj == nullptr)
		return;

	for (auto component : sobj->GetAllComponents())
		visitor->Visit(component);

	AddController(sobj);
}

void Attribute::AddController(Ptr<SObj> sobj) {
	auto item = new QWidget;
	tbox->insertItem(tbox->count(), item, "Controller");

	auto grid = Grid::New(item);
	item2grid[item] = grid;

	grid->AddTitle("[ Component ]");

	const int componentNum = 5;
	vector<string> componentNames{ "Camera" , "Geometry", "Light", "Material", "Transform" };

	vector<function<Ptr<Component>()>> componentGenFuncs{
		[=]()->Ptr<Component> { return CmptCamera::New(nullptr); },
		[=]()->Ptr<Component> { return CmptGeometry::New(nullptr, nullptr); },
		[=]()->Ptr<Component> { return CmptLight::New(nullptr, nullptr); },
		[=]()->Ptr<Component> { return CmptMaterial::New(nullptr, nullptr); },
		[=]()->Ptr<Component> { return CmptTransform::New(nullptr); },
	};

	vector<function<Ptr<Component>()>> componentDelFuncs{
		[=]()->Ptr<Component> {
			auto component = sobj->GetComponent<CmptCamera>();
			sobj->DetachComponent<CmptCamera>();
			return component;
		},
		[=]()->Ptr<Component> {
			auto component = sobj->GetComponent<CmptGeometry>();
			sobj->DetachComponent<CmptGeometry>();
			return component;
		},
		[=]()->Ptr<Component> {
			auto component = sobj->GetComponent<CmptLight>();
			sobj->DetachComponent<CmptLight>();
			return component;
		},
		[=]()->Ptr<Component> {
			auto component = sobj->GetComponent<CmptMaterial>();
			sobj->DetachComponent<CmptMaterial>();
			return component;
		},
		[=]()->Ptr<Component> {
			auto component = sobj->GetComponent<CmptTransform>();
			sobj->DetachComponent<CmptTransform>();
			return component;
		},
	};

	map<std::string, function<Ptr<Component>()> > componentGenMap;
	map<std::string, function<Ptr<Component>()> > componentDelMap;
	for (int i = 0; i < componentNum; i++) {
		componentGenMap[componentNames[i]] = componentGenFuncs[i];
		componentDelMap[componentNames[i]] = componentDelFuncs[i];
	}

	auto addBtnSlot = [=](const string& item) {
		auto target = componentGenMap.find(item);
		if (target == componentGenMap.cend())
			return;

		auto component = target->second();
		if (sobj->HaveComponentSameTypeWith(component))
			return;

		sobj->AttachComponent(component);
		visitor->Visit(component);
	};

	auto delBtnSlot = [=](const string& item) {
		auto target = componentDelMap.find(item);
		if (target == componentDelMap.cend())
			return;

		auto component = target->second();
		if (component == nullptr)
			return;

		auto target2 = componentType2item.find(typeid(*component));
		if (target2 == componentType2item.cend())
			return;

		auto tboxItem = target2->second;

		componentType2item.erase(target2);

		item2grid.erase(tboxItem);

		tbox->removeItem(tbox->indexOf(tboxItem));
		delete tboxItem;
	};

	grid->AddComboBox("- Add Component", componentNames, "Add", addBtnSlot);
	grid->AddComboBox("- Delete Component", componentNames, "Delete", delBtnSlot);
}
