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
	class Kinematic : public PlanarEmbedding {
	public:
		Kinematic(Ptr<TriMesh> triMesh) : PlanarEmbedding(triMesh) { }
		Kinematic(Ptr<Parameterization> embd) : PlanarEmbedding(embd) { }
		Kinematic(Ptr<Preprocessing> preproc) : PlanarEmbedding(preproc) { }
	public:
		static const Ptr<Kinematic> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Kinematic>(triMesh);
		}
		static const Ptr<Kinematic> New(Ptr<Parameterization> embd) {
			return Ubpa::New<Kinematic>(embd);
		}
		static const Ptr<Kinematic> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<Kinematic>(preproc);
		}
	public:
		bool Run();
		bool Iterate();

	protected:
		void InitPara();
		void InitSolver();
		void UpdateTriMesh();

	private:
		void calc_accel();

		void Tutte();
		void Pre_calculate();

		void Energysource();

		void local_coordinate_inverse(int t, Matrix2d& p);

		double delta_t, m;

		VectorXd area_weight;

		std::vector<Matrix2d> source_p;

		std::vector<int> F0, F1, F2;

		std::vector<double> energy_area_process;
		double energy_pre, energy_cur;
		double energy_uniform, energy_area;

		VectorXd input_mesh_pos;
		VectorXd position, velocity, acceleration;

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

