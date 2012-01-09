#include "OpenVolumeMeshBaseKernel.hh"

void OpenVolumeMeshBaseKernel::property_stats(std::ostream& _ostr) const {
	const OpenVolumeMeshPropertyContainer::Properties& vps =
			vprops_.properties();
	const OpenVolumeMeshPropertyContainer::Properties& hps =
			heprops_.properties();
	const OpenVolumeMeshPropertyContainer::Properties& eps =
			eprops_.properties();
	const OpenVolumeMeshPropertyContainer::Properties& fps =
			fprops_.properties();
	const OpenVolumeMeshPropertyContainer::Properties& hfps =
			hfprops_.properties();
	const OpenVolumeMeshPropertyContainer::Properties& cps =
			cprops_.properties();
	const OpenVolumeMeshPropertyContainer::Properties& mps =
			mprops_.properties();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;

	_ostr << vprops_.size() << " vprops:\n";
	for (it = vps.begin(); it != vps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
	_ostr << heprops_.size() << " heprops:\n";
	for (it = hps.begin(); it != hps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
	_ostr << eprops_.size() << " eprops:\n";
	for (it = eps.begin(); it != eps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
	_ostr << fprops_.size() << " fprops:\n";
	for (it = fps.begin(); it != fps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
	_ostr << hfprops_.size() << " hfprops:\n";
	for (it = hfps.begin(); it != hfps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
	_ostr << cprops_.size() << " cprops:\n";
	for (it = cps.begin(); it != cps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
	_ostr << mprops_.size() << " mprops:\n";
	for (it = mps.begin(); it != mps.end(); ++it) {
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
	}
}

void OpenVolumeMeshBaseKernel::vprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& vps =
			vprops_.properties();
	for (it = vps.begin(); it != vps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}
}

void OpenVolumeMeshBaseKernel::heprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& hps =
			heprops_.properties();
	for (it = hps.begin(); it != hps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}
}

void OpenVolumeMeshBaseKernel::eprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& eps =
			eprops_.properties();
	for (it = eps.begin(); it != eps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}
}

void OpenVolumeMeshBaseKernel::fprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& fps =
			fprops_.properties();
	for (it = fps.begin(); it != fps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}
}

void OpenVolumeMeshBaseKernel::hfprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& hfps =
			hfprops_.properties();
	for (it = hfps.begin(); it != hfps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}
}

void OpenVolumeMeshBaseKernel::cprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& cps =
			cprops_.properties();
	for (it = cps.begin(); it != cps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}
}

void OpenVolumeMeshBaseKernel::mprop_stats(std::string& _string) const {
	_string.clear();

	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& mps =
			mprops_.properties();
	for (it = mps.begin(); it != mps.end(); ++it)
		if (*it == NULL)
			_string += "[deleted] \n";
		else {
			_string += (*it)->name();
			_string += "\n";
		}

}

void OpenVolumeMeshBaseKernel::vprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& vps =
			vprops_.properties();
	for (it = vps.begin(); it != vps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}

void OpenVolumeMeshBaseKernel::heprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& hps =
			heprops_.properties();
	for (it = hps.begin(); it != hps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}

void OpenVolumeMeshBaseKernel::eprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& eps =
			eprops_.properties();
	for (it = eps.begin(); it != eps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}

void OpenVolumeMeshBaseKernel::fprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& fps =
			fprops_.properties();
	for (it = fps.begin(); it != fps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}

void OpenVolumeMeshBaseKernel::hfprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& hfps =
			hfprops_.properties();
	for (it = hfps.begin(); it != hfps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}

void OpenVolumeMeshBaseKernel::cprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& cps =
			cprops_.properties();
	for (it = cps.begin(); it != cps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}

void OpenVolumeMeshBaseKernel::mprop_stats(std::ostream& _ostr) const {
	OpenVolumeMeshPropertyContainer::Properties::const_iterator it;
	const OpenVolumeMeshPropertyContainer::Properties& mps =
			mprops_.properties();
	for (it = mps.begin(); it != mps.end(); ++it)
		*it == NULL ? (void) (_ostr << "[deleted]" << "\n") : (*it)->stats(
				_ostr);
}
