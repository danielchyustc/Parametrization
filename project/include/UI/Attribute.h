#pragma once

#include <Basic/TypeMap.h>

#include <Basic/Ptr.h>

#include <qobject.h>
#include <qtoolbox.h>
#include <Engine/MeshEdit/Parameterization/ParaEnum.h>
#include <Engine/MeshEdit/Parameterization/Parameterization.h>

#include <map>

namespace Ubpa {
	class RawAPI_OGLW;
	class SObj;
	class Component;
	class Parameterization;
	class Grid;

	enum WeightType;
	enum BoundShape;
	enum ViewMode;
	enum PlanarMethod;
	enum SphericalMethod;

	class Attribute final {
	protected:
		Attribute();

	public:
		static Attribute* GetInstance() {
			static Attribute instance;
			return &instance;
		}

	public:
		void Init(QToolBox* tbox, RawAPI_OGLW* pOGLW);
		void SetSObj(Ptr<SObj> sobj);
		const Ptr<SObj> GetCurSObj() const { return curSObj.lock(); }
		void SetWeightType(WeightType weight_type) { weight_type_ = weight_type; }
		void SetBoundShape(BoundShape bound_shape) { bound_shape_ = bound_shape; }
		void SetViewMode(ViewMode view_mode) {
			view_mode_ = view_mode;
		}
		template<typename T, typename = std::enable_if_t<std::is_base_of_v<Component, T>>>
		void SetCurCmpt() {
			auto target = componentType2item.find(typeid(T));
			if (target == componentType2item.end())
				return;

			tbox->setCurrentWidget(target->second);
		}

	private:
		void AddController(Ptr<SObj> sobj);

	private:
		class ComponentVisitor;
		friend class ComponentVisitor;

		QToolBox* tbox;

		TypeMap<QWidget*> componentType2item;
		std::map<QWidget*, Ptr<Grid>> item2grid;

		Ptr<ComponentVisitor> visitor;

		Ptr<Parameterization> param_;
		PlanarMethod plane_method_;
		SphericalMethod sphere_method_;
		WPtr<SObj> curSObj;

		RawAPI_OGLW* pOGLW;

		WeightType weight_type_ = kEqual;
		BoundShape bound_shape_ = kCircle;
		ViewMode view_mode_ = kMesh;
	};
}
