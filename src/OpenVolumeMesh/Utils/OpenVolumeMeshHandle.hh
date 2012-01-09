/*
 * HandleType.hh
 *
 *  Created on: 27.06.2011
 *      Author: mike
 */

#ifndef OPENVOLUMEMESHHANDLE_HH_
#define OPENVOLUMEMESHHANDLE_HH_

// Define handle types in order to distinguish different entities by their indices
class OpenVolumeMeshHandle {
public:
    // Default constructor
	explicit OpenVolumeMeshHandle(int _idx) : idx_(_idx) {};

	OpenVolumeMeshHandle& operator=(int _idx) {
		idx_ = _idx;
		return *this;
	}

	OpenVolumeMeshHandle& operator=(const OpenVolumeMeshHandle& _idx) {
		idx_ = _idx.idx_;
		return *this;
	}

	inline bool is_valid() const { return idx_ != -1; }

	inline bool operator<(const OpenVolumeMeshHandle& _idx) const { return (this->idx_ < _idx.idx_); }

	inline bool operator<(int _idx) const { return idx_ < _idx; }

	inline int idx() const { return idx_; }

	void idx(const int& _idx) { idx_ = _idx; }

	void reset() { idx_ = -1; }

	operator int() const { return idx_; }

private:
	int idx_;
};

// Default entity handles

class VertexHandle   : public OpenVolumeMeshHandle { public: VertexHandle(int _idx = -1)   : OpenVolumeMeshHandle(_idx) {} };
class EdgeHandle     : public OpenVolumeMeshHandle { public: EdgeHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class FaceHandle     : public OpenVolumeMeshHandle { public: FaceHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class CellHandle     : public OpenVolumeMeshHandle { public: CellHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class HalfEdgeHandle : public OpenVolumeMeshHandle { public: HalfEdgeHandle(int _idx = -1) : OpenVolumeMeshHandle(_idx) {} };
class HalfFaceHandle : public OpenVolumeMeshHandle { public: HalfFaceHandle(int _idx = -1) : OpenVolumeMeshHandle(_idx) {} };

#endif /* OPENVOLUMEMESHHANDLE_HH_ */
