/*
 * Status.hh
 *
 *  Created on: Jun 10, 2011
 *      Author: kremer
 */

#ifndef STATUS_HH_
#define STATUS_HH_

class OpenVolumeMeshStatus {
public:

    // Default constructor
    OpenVolumeMeshStatus() : selected_(false), tagged_(false), deleted_(false) {}

    bool selected() const { return selected_; }

    bool tagged() const { return tagged_; }

    bool deleted() const { return deleted_; }

    void set_selected(bool _selected) { selected_ = _selected; }

    void set_tagged(bool _tagged) { tagged_ = _tagged; }

    void set_deleted(bool _deleted) { deleted_ = _deleted; }

private:

    bool selected_;

    bool tagged_;

    bool deleted_;
};

#endif /* STATUS_HH_ */
