#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <Engine/MeshEdit/Parameterization/ParaEnum.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class Preprocessing;
	enum ViewMode;

	// mesh boundary == 1
	class Parameterization : public HeapObj {
	public:
		Parameterization(Ptr<TriMesh> triMesh);
		Parameterization(Ptr<Parameterization> para);
		Parameterization(Ptr<Preprocessing> preproc);
	public:
		static const Ptr<Parameterization> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Parameterization>(triMesh);
		}
		static const Ptr<Parameterization> New(Ptr<Parameterization> para) {
			return Ubpa::New<Parameterization>(para);
		}
	public:
		virtual bool Init(Ptr<TriMesh> triMesh) { return true; }
		virtual void Clear() { return; }
		virtual bool Run() { return true; }
		virtual void UpdateTriMesh() { }
		virtual bool Iterate() { return true; }
		virtual bool IterateTillDone() { return true; }
		virtual void SetViewMode(ViewMode view_mode) { 
			view_mode_ = view_mode; 
			if(init_status_)
				UpdateTriMesh(); 
		}
		virtual void MoveBarycenterToOrigin() { }
		virtual void AreaCorrection() { }
		ViewMode GetViewMode() { return view_mode_; }
		Ptr<TriMesh> GetTriMesh() { return triMesh; }
		MatrixXd ParaMap() { return para_map_; }
		bool FinishStatus() { return finish_status_; }

	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P : public TPolygon<V, E, P> { };
	protected:
		bool InitHeMesh();
		Vector3d NormalAt(int i);
		double Cotan(V* vi, V* vj);
		double CotanWeight(int i, V* vi, int j, V* vj);
		double CotanWeight(int i, V* vi, V* vj);
		double CotanWeight(V* vi, V* vj);
		double CotanWeight(int i, int j);
		double SurfaceElement(V* vi);
		double SurfaceElement(int i) { return SurfaceElement(heMesh->Vertices().at(i)); }
		double Barycenter(RowVectorXd f);
		VectorXd Barycenter(MatrixXd f);
		double GaussCurvature(int i);
		double MeanCurvature(int i);
		double TriArea(int t);
		vecf3 TriNormal(int t);
		double VoronoiArea(V* vi);
		double VoronoiArea(int i) { return VoronoiArea(heMesh->Vertices().at(i)); }
		double TotalArea();
		vecf3 CircumCenter(V* vi, V* vj);
		void BuildCotanWeightMatrix();
	protected:
		friend class Minimize;

		Ptr<TriMesh>			triMesh;
		Ptr<HEMesh<V>>			heMesh;
		Ptr<Preprocessing>		preproc;
		size_t					nV;
		size_t					nT;
		ViewMode				view_mode_;
		MatrixXd				para_map_;
		SparseMatrix<double, RowMajor>	CotW;
		int						iter_count_ = 0;
		bool					finish_status_ = false;
		bool					init_status_ = false;
		bool					cot_mat_available_ = false;
	};
}
