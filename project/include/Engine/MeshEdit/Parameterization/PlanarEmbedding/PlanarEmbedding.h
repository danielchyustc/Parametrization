#pragma once

#include <Engine/MeshEdit/Parameterization/Parameterization.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class Preprocessing;
	class Parameterization;

	// mesh boundary == 1
	class PlanarEmbedding : public Parameterization {
	public:
		PlanarEmbedding(Ptr<TriMesh> triMesh);
		PlanarEmbedding(Ptr<Parameterization> embd);
		PlanarEmbedding(Ptr<Preprocessing> preproc);
	public:
		static const Ptr<PlanarEmbedding> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<PlanarEmbedding>(triMesh);
		}
		static const Ptr<PlanarEmbedding> New(Ptr<Parameterization> embd) {
			return Ubpa::New<PlanarEmbedding>(embd);
		}
		static const Ptr<PlanarEmbedding> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<PlanarEmbedding>(preproc);
		}
	public:
		virtual void Clear() { return; }
		virtual bool Init(Ptr<TriMesh> triMesh) { return true; }
		virtual bool Run();
		virtual bool Iterate() { return true; }
		virtual bool IterateTillDone() { return true; }
		void SetViewMode(ViewMode view_mode) {
			view_mode_ = view_mode;
			if (init_status_)
				UpdateTriMesh();
		}

	protected:
		void UpdateParaMap() { }
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
	};
}
