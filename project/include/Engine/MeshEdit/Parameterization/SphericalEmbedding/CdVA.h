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
	class CdVA : public SphericalEmbedding {
	public:
		CdVA(Ptr<TriMesh> triMesh) : SphericalEmbedding(triMesh) { }
		CdVA(Ptr<Parameterization> embd) : SphericalEmbedding(embd) { }
		CdVA(Ptr<Preprocessing> preproc) : SphericalEmbedding(preproc) { }
	public:
		static const Ptr<CdVA> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<CdVA>(triMesh);
		}
		static const Ptr<CdVA> New(Ptr<Parameterization> embd) {
			return Ubpa::New<CdVA>(embd);
		}
		static const Ptr<CdVA> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<CdVA>(preproc);
		}
	public:
		bool Run();
		bool Iterate();
		bool IterateTillDone();

	private:
		void InitSolver();
		void UpdateParaMap();
		double Energy();

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