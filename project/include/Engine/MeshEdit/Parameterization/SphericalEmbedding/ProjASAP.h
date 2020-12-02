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
	class ProjASAP : public SphericalEmbedding {
	public:
		ProjASAP(Ptr<TriMesh> triMesh) : SphericalEmbedding(triMesh) { }
		ProjASAP(Ptr<Parameterization> embd) : SphericalEmbedding(embd) { }
		ProjASAP(Ptr<Preprocessing> preproc) : SphericalEmbedding(preproc) { }
	public:
		static const Ptr<ProjASAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<ProjASAP>(triMesh);
		}
		static const Ptr<ProjASAP> New(Ptr<Parameterization> embd) {
			return Ubpa::New<ProjASAP>(embd);
		}static const Ptr<ProjASAP> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<ProjASAP>(preproc);
		}
	public:
		bool Run();
		bool Iterate() { return true; }
		bool IterateTillDone() { return true; }

	private:
		void InitParaMap();
		void UpdateParaMap() { return; }

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
		Matrix<double, 3, 2>* flat_tri_;
	};
}