#pragma once
// Composed Hessian Metric
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
	class CHM : public PlanarEmbedding {
	public:
		CHM(Ptr<TriMesh> triMesh) : PlanarEmbedding(triMesh) { }
		CHM(Ptr<Parameterization> embd) : PlanarEmbedding(embd) { }
		CHM(Ptr<Preprocessing> preproc) : PlanarEmbedding(preproc) { }
	public:
		static const Ptr<CHM> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<CHM>(triMesh);
		}
		static const Ptr<CHM> New(Ptr<Parameterization> embd) {
			return Ubpa::New<CHM>(embd);
		}
		static const Ptr<CHM> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<CHM>(preproc);
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
		void sym_dir_energy();
		void sym_dir_energy(VectorXd pos, double& e);
		void negative_gradient();

		void solve_velocity();
		void max_step(const VectorXd& xx, const VectorXd& dd, double& step);
		void backtracking_line_search(const VectorXd& x, const VectorXd& d, const VectorXd& negetive_grad, double& alpha);

		MatrixXd source_lc;
		double energy;
		VectorXd neg_gradient;

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

