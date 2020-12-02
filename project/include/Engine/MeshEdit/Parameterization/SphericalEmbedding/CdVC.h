#pragma once

#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SphericalEmbedding.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class Preprocessing;
	enum ViewMode;

	// mesh boundary == 0
	// spherical barrycentric embedding
	class CdVC : public SphericalEmbedding {
	public:
		CdVC(Ptr<TriMesh> triMesh) : SphericalEmbedding(triMesh) { }
		CdVC(Ptr<Parameterization> embd) : SphericalEmbedding(embd) { }
		CdVC(Ptr<Preprocessing> preproc) : SphericalEmbedding(preproc) { }
	public:
		static const Ptr<CdVC> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<CdVC>(triMesh);
		}
		static const Ptr<CdVC> New(Ptr<Parameterization> embd) {
			return Ubpa::New<CdVC>(embd);
		}
		static const Ptr<CdVC> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<CdVC>(preproc);
		}
	public:
		bool Run();
		bool Iterate();
		bool IterateTillDone();

	private:
		void InitSolver();
		void UpdateParaMap();
		double Energy(RowVectorXd f);
		double Energy(MatrixXd f);

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
	private:
		double					energy_;
		double					epsilon = 0.001;
		SparseMatrix<double, RowMajor>	L;
		SparseLU<SparseMatrix<double>> solver;
	};
}