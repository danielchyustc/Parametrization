#pragma once

#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/PlanarEmbedding.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <UHEMesh/HEMesh.h>
#include <UI/Attribute.h>
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
	class ARAP : public PlanarEmbedding {
	public:
		ARAP(Ptr<TriMesh> triMesh) : PlanarEmbedding(triMesh) { }
		ARAP(Ptr<Parameterization> embd) : PlanarEmbedding(embd) { }
		ARAP(Ptr<Preprocessing> preproc) : PlanarEmbedding(preproc) { }
	public:
		static const Ptr<ARAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<ARAP>(triMesh);
		}
		static const Ptr<ARAP> New(Ptr<Parameterization> embd) {
			return Ubpa::New<ARAP>(embd);
		}
		static const Ptr<ARAP> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<ARAP>(preproc);
		}
	public:
		bool Run();
		bool Iterate();

	protected:
		void InitPara() { }
		void InitFlatTri();
		void InitParaSolver();
		void UpdateParaMap();
		void UpdateTriMesh();

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
		Matrix<double, 3, 2>				*flat_tri_;
		SparseLU<SparseMatrix<double>>	para_solver_;
	};
}

