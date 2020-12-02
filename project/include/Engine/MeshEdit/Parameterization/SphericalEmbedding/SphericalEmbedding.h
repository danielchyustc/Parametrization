#pragma once

#include <Engine/MeshEdit/Parameterization/Parameterization.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class Parameterization;
	class Preprocessing;
	enum ViewMode;

	// mesh boundary == 0
	// spherical embedding
	class SphericalEmbedding : public Parameterization {
	public:
		SphericalEmbedding(Ptr<TriMesh> triMesh);
		SphericalEmbedding(Ptr<Parameterization> embd);
		SphericalEmbedding(Ptr<Preprocessing> preproc);
	public:
		static const Ptr<SphericalEmbedding> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<SphericalEmbedding>(triMesh);
		}
		static const Ptr<SphericalEmbedding> New(Ptr<Parameterization> embd) {
			return Ubpa::New<SphericalEmbedding>(embd);
		}
		static const Ptr<SphericalEmbedding> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<SphericalEmbedding>(preproc);
		}
	public:
		virtual void Clear() { return; }
		virtual bool Init(Ptr<TriMesh> triMesh) { return true; }
		virtual bool Run();
		virtual bool Iterate();
		virtual bool IterateTillDone();
		void UpdateTriMesh();
		void SetViewMode(ViewMode view_mode) {
			view_mode_ = view_mode;
			if (init_status_)
				UpdateTriMesh();
		}
		void MoveBarycenterToOrigin();
		void AreaCorrection();

	protected:
		void InitParaMap();
		virtual void UpdateParaMap() { }
		Matrix3d Projection(Vector3d v);
		Vector2d StereoProjection(Vector3d x) {
			return Vector2d(x[0] / (1 + x[2]), x[1] / (1 + x[2]));
		}

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
		double radius = 1;
	};
}