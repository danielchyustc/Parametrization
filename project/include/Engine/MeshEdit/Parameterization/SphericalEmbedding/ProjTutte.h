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
	class ProjTutte : public SphericalEmbedding {
	public:
		ProjTutte(Ptr<TriMesh> triMesh) : SphericalEmbedding(triMesh) { }
		ProjTutte(Ptr<Parameterization> embd) : SphericalEmbedding(embd) { }
		ProjTutte(Ptr<Preprocessing> preproc) : SphericalEmbedding(preproc) { }
	public:
		static const Ptr<ProjTutte> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<ProjTutte>(triMesh);
		}
		static const Ptr<ProjTutte> New(Ptr<Parameterization> embd) {
			return Ubpa::New<ProjTutte>(embd);
		}static const Ptr<ProjTutte> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<ProjTutte>(preproc);
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
	};
}