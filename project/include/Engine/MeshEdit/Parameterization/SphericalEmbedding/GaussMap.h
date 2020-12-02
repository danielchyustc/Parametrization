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
	class GaussMap : public SphericalEmbedding {
	public:
		GaussMap(Ptr<TriMesh> triMesh) : SphericalEmbedding(triMesh) { }
		GaussMap(Ptr<Parameterization> embd) : SphericalEmbedding(embd) { }
		GaussMap(Ptr<Preprocessing> preproc) : SphericalEmbedding(preproc) { }
	public:
		static const Ptr<GaussMap> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<GaussMap>(triMesh);
		}
		static const Ptr<GaussMap> New(Ptr<Parameterization> embd) {
			return Ubpa::New<GaussMap>(embd);
		}static const Ptr<GaussMap> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<GaussMap>(preproc);
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