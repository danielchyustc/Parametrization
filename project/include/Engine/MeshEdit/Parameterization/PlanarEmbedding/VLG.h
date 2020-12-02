#pragma once

#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/PlanarEmbedding.h>
#include <Engine/MeshEdit/Parameterization/PardisoSolver.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class Preprocessing;
	class Paramaterize;
	enum ViewMode;
	// mesh boundary == 1
	class VLG : public PlanarEmbedding {
	public:
		VLG(Ptr<TriMesh> triMesh) : PlanarEmbedding(triMesh) { }
		VLG(Ptr<Parameterization> embd) : PlanarEmbedding(embd) { }
		VLG(Ptr<Preprocessing> preproc) : PlanarEmbedding(preproc) { }
	public:
		static const Ptr<VLG> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<VLG>(triMesh);
		}
		static const Ptr<VLG> New(Ptr<Parameterization> embd) {
			return Ubpa::New<VLG>(embd);
		}
		static const Ptr<VLG> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<VLG>(preproc);
		}
	public:
		bool Run();
		bool Iterate();

	protected:
		void InitPara();
		void UpdateTriMesh();

	private:
		void Tutte();

		void source_local_coordinate();
		void target_local_coordinate();

		void solve_velocity();

		MatrixXd source_lc;
		MatrixXd target_lc;

		std::vector<int> F0, F1, F2;

		VectorXd area_weight;

		VectorXd input_mesh_pos;
		VectorXd position, velocity;

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
	};
}

