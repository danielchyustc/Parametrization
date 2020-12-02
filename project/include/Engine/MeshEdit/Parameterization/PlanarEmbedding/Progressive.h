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
	class Progressive : public PlanarEmbedding {
	public:
		Progressive(Ptr<TriMesh> triMesh) : PlanarEmbedding(triMesh) { }
		Progressive(Ptr<Parameterization> embd) : PlanarEmbedding(embd) { }
		Progressive(Ptr<Preprocessing> preproc) : PlanarEmbedding(preproc) { }
	public:
		static const Ptr<Progressive> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Progressive>(triMesh);
		}
		static const Ptr<Progressive> New(Ptr<Parameterization> embd) {
			return Ubpa::New<Progressive>(embd);
		}
		static const Ptr<Progressive> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<Progressive>(preproc);
		}
	public:
		bool Run();
		bool Iterate();

	protected:
		void InitPara();
		void InitSolver();
		void UpdateTriMesh();

	private:
		void calc_gradient_norm(const VectorXd& x);

		void recover_to_src();

		void Update_source_same_t();

		void Tutte();
		void Pre_calculate();

		void CM();
		void SLIM();

		void Energy(const VectorXd& x, double& energy);
		void Energysource();

		double newton_equation(const double& a, const double& b, const double& K);

		void backtracking_line_search(const VectorXd& x, const VectorXd& d, const VectorXd& negetive_grad, double& alpha);

		void local_coordinate_inverse(int t, Matrix2d& p);

		void max_step(const VectorXd& xx, const VectorXd& dd, double& step);

		double Intp_T_Min;
		double changetocm_flag;
		bool flag_1, flag_2;
		double convgence_con_rate;
		double time_consumption;
		int MAX_ITER_NUM;
		double originmesh_area_sqrt;

		VectorXd negative_grad_norm;
		double g_norm;

		VectorXd area, area_uniform, area_src;

		std::vector<Matrix2d> source_p, update_p;

		std::vector<int> F0, F1, F2;

		PardisoSolver* pardiso;
		std::vector<int> pardiso_ia, pardiso_i, pardiso_ja;
		std::vector<double> pardiso_a, pardiso_b;

		std::vector<double> energy_area_process;
		double energy_pre, energy_cur;
		double energy_uniform, energy_area;

		double bound_distortion_K;
		VectorXd input_mesh_pos;
		VectorXd position_of_mesh;

		int slim_iter_num, cm_iter_num, sum_iter_num;

		double conv_percent;

		std::vector<int> id_h00, id_h01, id_h02, id_h03, id_h04, id_h05;
		std::vector<int> id_h11, id_h12, id_h13, id_h14, id_h15;
		std::vector<int> id_h22, id_h23, id_h24, id_h25;
		std::vector<int> id_h33, id_h34, id_h35;
		std::vector<int> id_h44, id_h45;
		std::vector<int> id_h55;

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

