#pragma once

namespace Ubpa {
	enum WeightType {
		kEqual,
		kCotangent,
	};

	enum BoundShape {
		kCircle,
		kSquare,
	};

	enum ViewMode {
		kMesh,
		kSphere,
		kPlane
	};

	enum PlanarMethod {
		kTutteEmbedding,
		kASAP,
		kARAP,
		kProg,
		kKine,
		kVLG,
		kCHM
	};

	enum SphericalMethod {
		kGaussMap,
		kProjTutte,
		kProjASAP,
		kSBE,
		kSCE,
		kCdVA,
		kCdVB,
		kCdVC,
		kCdVARAP,
		kSARAP
	};
}
