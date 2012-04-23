#include "ofxProCamSolver/CorrespondenceSet.h"

namespace ofxProCamSolver {

	//---------
	template<typename T>
	void CorrespondenceSet_<T>::save(string filename) {
		if (filename=="")
			filename = ofSystemSaveDialog("correspondences", "Save ofxProCamSolver::CorrespondenceSet").getPath();
		if (filename=="") {
			ofLogError("ofxProCamSolver") << "No file selected for save, aborting";
			return;
		}

		//count
		uint32_t count = this->size();
		ofstream file;
		file.open(filename, ios::binary);
		file.write((char*)&count, sizeof(uint32_t));

		//items
		CorrespondenceSet_<T>::iterator it;
		for (it = this->begin(); it != this->end(); it++)
			file.write((char*)&*it, sizeof(Correspondence_<T>));
		file.close();
	}

	
	//---------
	template<typename T>
	void CorrespondenceSet_<T>::load(string filename) {
		if (filename=="")
			filename = ofSystemLoadDialog("Load ofxProCamSolver::CorrespondenceSet").getPath();
		if (filename=="") {
			ofLogError("ofxProCamSolver") << "No file selected for load, aborting";
			return;
		}
		
		//count
		uint32_t count;
		ifstream file;
		file.open(filename, ios::binary);
		file.read((char*)&count, sizeof(uint32_t));
		this->resize(count);

		//items
		CorrespondenceSet_<T>::iterator it;
		for (it = this->begin(); it != this->end(); it++)
			file.read((char*)&(*it), sizeof(Correspondence_<T>));
	}

	//---------
	template<typename T>
	set<int> CorrespondenceSet_<T>::getViewIndices() const {
		set<int> indices;
		CorrespondenceSet_<T>::const_iterator it;
		for (it = this->begin(); it != this->end(); it++) {
			if (indices.count(it->cameraID1) == 0)
				indices.insert(it->cameraID1);
			if (indices.count(it->cameraID2) == 0)
				indices.insert(it->cameraID2);
		}

		return indices;
	}

#ifdef HAVE_OFXGRAYCODE
	using namespace ofxGraycode;

	//---------
	template<typename T>
	void CorrespondenceSet_<T>::add(const ofxGraycode::DataSet &dataSet, int cameraIndex, int projectorIndex) {
		if (cameraIndex==projectorIndex) {
			ofLogError("ofxProCamSolver") << "addDataSet: cameraIndex == projectorIndex. All camera / projector indices must be unique. e.g. a scene with 3 cameras and 2 projectors would have indices 0,1,2,3,4";
			return;
		}

		//build a set of camera finds for each projector point
		//output average, sdev finds per projector point
		map<uint32_t, ProjectorPixel> ppixels;
		const uint32_t *data = dataSet.getData().getPixels();
		const uint32_t *distance = dataSet.getDistance().getPixels();
		const uint8_t* active = dataSet.getActive().getPixels();
		for (uint32_t i=0; i<dataSet.size(); i++, data++, distance++, active++) {
			if (*active) {
				if (ppixels.count(*data) != 0)
					ppixels.insert(pair<uint32_t, ProjectorPixel>(*data, ProjectorPixel(i, *data, *distance)));
				else
					ppixels[*data].addCameraFind(i, *distance);
			}
		}
	
		//create new ofxProCamSolver correspondence set
		this->reserve(this->size() + ppixels.size());

		//loop through projector pixels and add weighted average to correspondence set
		//save correspondence set
		map<uint32_t, ProjectorPixel>::const_iterator it;
		for (it = ppixels.begin(); it != ppixels.end(); it++)
			this->push_back(Correspondence_<T>(cameraIndex, it->second.projector, projectorIndex, it->second.getCameraMean()));
	}
#endif

	template class CorrespondenceSet_<double>;
	template class CorrespondenceSet_<float>;
}