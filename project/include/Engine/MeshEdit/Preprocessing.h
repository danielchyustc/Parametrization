#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <UI/Attribute.h>

namespace Ubpa {
	class TriMesh;
	class Paramaterize;

	class Preprocessing : public HeapObj {
	public:
		Preprocessing(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<Preprocessing> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Preprocessing>(triMesh);
		}
	public:
		// clear cache data
		void Clear();

		// init cache data (eg. half-edge structure) for Run()
		bool Init(Ptr<TriMesh> triMesh);

		size_t NumVertices() { return nV; }
		size_t NumEdges() { return nE; }
		size_t NumTriangles() { return nT; }
		size_t NumBoundaries() { return nB; }
		size_t EulerChar() { return euler_char; }
		size_t Genus() { return genus; }
		bool IsTriMesh() { return isTriMesh; }
		bool IsClosed() { return isClosed; }

	private:
		bool Glue(Ptr<TriMesh> triMesh);

	private:
		size_t	nV;
		size_t	nE;
		size_t	nT;
		size_t	nB;
		size_t	euler_char;
		size_t	genus;
		bool	isTriMesh;
		bool	isClosed;

	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P :public TPolygon<V, E, P> { };
	public:
		Ptr<TriMesh> GetTriMesh() { return triMesh; }
		Ptr<HEMesh<V>> GetHEMesh() { return heMesh; }
	private:
		friend class Paramaterize;

		Ptr<TriMesh> triMesh;
		Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh
	};
}
