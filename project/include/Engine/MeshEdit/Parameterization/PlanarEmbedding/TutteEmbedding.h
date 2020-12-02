#pragma once

#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/PlanarEmbedding.h>
#include <Engine/MeshEdit/Preprocessing.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <UI/Attribute.h>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class Preprocessing;

	// mesh boundary == 1
	class TutteEmbedding : public PlanarEmbedding {
	public:
		TutteEmbedding(Ptr<TriMesh> triMesh) : PlanarEmbedding(triMesh) { }
		TutteEmbedding(Ptr<Parameterization> embd) : PlanarEmbedding(embd) { }
		TutteEmbedding(Ptr<Preprocessing> preproc) : PlanarEmbedding(preproc) { }
	public:
		static const Ptr<TutteEmbedding> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<TutteEmbedding>(triMesh);
		}
		static const Ptr<TutteEmbedding> New(Ptr<Parameterization> embd) {
			return Ubpa::New<TutteEmbedding>(embd);
		}
		static const Ptr<TutteEmbedding> New(Ptr<Preprocessing> preproc) {
			return Ubpa::New<TutteEmbedding>(preproc);
		}
	public:
		bool Run();
		void SetWeightType(WeightType weight_type) { weight_type_ = weight_type; }
		void SetBoundShape(BoundShape bound_shape) { bound_shape_ = bound_shape; }

	private:
		void UpdateParaMap();

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
		WeightType weight_type_;
		BoundShape bound_shape_;
	};
}
