// $Id: $
//====================================================================
//  Based on ZPlanarTracker module from F. Gaede

//  Tracking detector to describe the FCC-ee IDEA vertex detector barrel.
//  The vertex detector is assembled of stave structures which feature 
//  support and readout (flex) elements. Each stave features multiple 
//  individual modules, that consist of sensitive and insensitive
//  sensor elements. From o1_v02 it is possible to define curved 
//  sensor and support elements. 
//--------------------------------------------------------------------
//
//  Author     : Armin Ilg
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"
#include <exception>
#include "DDRec/Surface.h"

using namespace std;

using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::Tube;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Volume;
using dd4hep::getAttrOrDefault;
using dd4hep::_toString;

using dd4hep::rec::volSurfaceList;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolPlane;


static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens)  {

    xml_det_t    x_det = e;
    int         m_id=0;
    std::string  det_name  = x_det.nameStr();

    DetElement sdet( det_name, x_det.id()  ) ;
    PlacedVolume pv;

    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

    envelope.setVisAttributes(theDetector, x_det.visStr());
    sens.setType("tracker");

    // Struct to support multiple readout or support layers (both using this struct)
    struct componentsStruct{
        string name;
        double r;
        double offset;
        double length;
        double z_offset;
        vector<double> thicknesses;
        vector<double> offsets; 
        vector<double> z_offsets; 
        vector<double> rs; 
        vector<Volume> volumes;
        vector<bool> isCurved;
    };

    // Struct to support end-of-stave structures
    struct endOfStaveStruct{
        string name;
        double r;
        double offset;
        double z_offset;
        vector<double> thicknesses;
        vector<double> lengths;
        vector<double> dzs;
        vector<double> offsets; 
        vector<double> z_offsets; 
        vector<double> rs; 
        vector<Volume> volumes;
        vector<int> side;
        vector<bool> isCurved;
    };

    // --- Module information struct ---
    struct stave_information{
        string name;
        double motherVolThickness;
        double motherVolWidth;

        vector<componentsStruct> components_vec;
        vector<endOfStaveStruct> endOfStaves_vec;
        double sensor_r;
        double sensor_offset;
        double sensor_thickness;
        Material sensor_material;
        vector<bool> sensor_sensitives;
        vector<double> sensor_xmin;
        vector<double> sensor_xmax;
        vector<double> sensor_ymin;
        vector<double> sensor_ymax;
        vector<double> sensor_phi_offsets;
        vector<bool> isCurved;
        double sensor_width;
        double sensor_length;
        double stave_dr;
        double stave_r;
        string staveVis;
        vector<Volume> sensor_volumes;
    };
    list<stave_information> stave_information_list;

    // --- Collect stave(s) information
    for(xml_coll_t mi(x_det,_U(stave)); mi; ++mi, ++m_id)  {
        xml_comp_t x_stave   = mi;

        stave_information m;
        m.name = x_stave.nameStr();
        m.staveVis = x_stave.visStr();

        m.stave_dr = x_stave.dr(0); // Offset for every second module in r
        m.stave_r = x_stave.r(0); // Radius of stave

        m.motherVolThickness = getAttrOrDefault(x_stave, _Unicode(motherVolThickness), double(0.0));
        m.motherVolWidth     = getAttrOrDefault(x_stave, _Unicode(motherVolWidth), double(0.0));

        // Components
        xml_coll_t c_components(x_stave,_U(components));
        for(c_components.reset(); c_components; ++c_components){
            componentsStruct components;
            components.name = xml_comp_t(c_components).nameStr();
            components.r = xml_comp_t(c_components).r(0);
            components.length = xml_comp_t(c_components).length();
            components.offset = xml_comp_t(c_components).offset(0);
            components.z_offset = xml_comp_t(c_components).z_offset(0);


            xml_coll_t c_component(c_components,_U(component));
            int iComponent = 0;
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                components.thicknesses.push_back(component.thickness());
                components.offsets.push_back(component.offset(0));
                components.z_offsets.push_back(component.z_offset(0));
                components.rs.push_back(component.r(0));

                bool isCurved = getAttrOrDefault(component, _Unicode(isCurved), bool(false));
                components.isCurved.push_back(isCurved);
                Volume ele_vol;

                if(isCurved){
                    double rmin = m.stave_r + components.r+components.rs.back();
                    double half_width = component.width()/(2.*M_PI*rmin)*(2.0*M_PI)/2.;
                    double phi_offset = getAttrOrDefault(component, _Unicode(phi_offset), double(0.0));
                    Tube ele_box = Tube(rmin, rmin+components.thicknesses.back(), components.length, -half_width + phi_offset, half_width + phi_offset);
                    ele_vol = Volume(components.name + _toString(iComponent, "_%d"), ele_box, theDetector.material(component.materialStr()));
                }
                else{
                    Box ele_box = Box(component.thickness()/2., component.width()/2., components.length);
                    ele_vol = Volume(components.name + _toString(iComponent, "_%d"), ele_box, theDetector.material(component.materialStr()));      
                }
                ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());  
                components.volumes.push_back(ele_vol);

                iComponent++;
            }
            m.components_vec.push_back(components);
        }

        // End of stave structures
        xml_coll_t c_endOfStave(x_stave,_U(end_z));
        for(c_endOfStave.reset(); c_endOfStave; ++c_endOfStave){
            endOfStaveStruct endOfStave;    
            endOfStave.offset = xml_comp_t(c_endOfStave).offset(0);
            endOfStave.z_offset = xml_comp_t(c_endOfStave).z_offset(0);
            endOfStave.r =      xml_comp_t(c_endOfStave).r(0);
            endOfStave.name =   xml_comp_t(c_endOfStave).nameStr();

            xml_coll_t c_component = xml_coll_t(c_endOfStave,_U(component));
            int iEndOfStave = 0;
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                endOfStave.thicknesses.push_back(component.thickness());
                endOfStave.dzs.push_back(component.dz(0));
                endOfStave.offsets.push_back(component.offset(0));
                endOfStave.z_offsets.push_back(component.z_offset(0));
                endOfStave.lengths.push_back(component.length());
                endOfStave.rs.push_back(component.r(0));
                
                int side = getAttrOrDefault(component, _Unicode(side), int(0));

                endOfStave.side.push_back(side); // 0 for both sides (default), +1 for +z side, -1 for -z side
                bool isCurved = getAttrOrDefault(component, _Unicode(isCurved), bool(false));
                endOfStave.isCurved.push_back(isCurved);
                Volume ele_vol;

                if(isCurved){
                    double rmin = m.stave_r + endOfStave.r+endOfStave.rs.back();
                    double half_width = component.width()/(2.*M_PI*rmin)*(2.0*M_PI)/2.;
                    double phi_offset = getAttrOrDefault(component, _Unicode(phi_offset), double(0.0));
                    Tube ele_box = Tube(rmin, rmin+endOfStave.thicknesses.back(), component.length()/2., -half_width + phi_offset, half_width + phi_offset);
                    ele_vol = Volume(endOfStave.name + _toString(iEndOfStave, "_%d"), ele_box, theDetector.material(component.materialStr()));
                }
                else{
                    Box ele_box = Box(component.thickness()/2., component.width()/2., component.length()/2.);
                    ele_vol = Volume(endOfStave.name + _toString(iEndOfStave, "_%d"), ele_box, theDetector.material(component.materialStr()));            
                }
                ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());

                endOfStave.volumes.push_back(ele_vol);
                iEndOfStave++;
            }
            m.endOfStaves_vec.push_back(endOfStave);
        }

        // Sensor
        xml_coll_t c_sensor(x_stave,_U(sensor));
        m.sensor_r = xml_comp_t(c_sensor).r(0);
        m.sensor_offset = xml_comp_t(c_sensor).offset(0);
        m.sensor_thickness = xml_comp_t(c_sensor).thickness();
        m.sensor_material = theDetector.material(xml_comp_t(c_sensor).materialStr());
        xml_coll_t c_component = xml_coll_t(c_sensor,_U(component));

        int iSensor = 0;
        for(c_component.reset(); c_component; ++c_component){
            xml_comp_t component = c_component;
            m.sensor_sensitives.push_back(component.isSensitive());
            m.sensor_xmin.push_back(component.xmin());
            m.sensor_xmax.push_back(component.xmax());
            m.sensor_ymin.push_back(component.ymin());
            m.sensor_ymax.push_back(component.ymax());

            m.sensor_phi_offsets.push_back(getAttrOrDefault(component, _Unicode(phi_offset), double(0.0)) );

            bool isCurved = getAttrOrDefault(component, _Unicode(isCurved), bool(false));
            m.isCurved.push_back(isCurved);

            // Already create volumes for all sensor components as this is independent of number of sensors per layer
            Volume ele_vol;
            if(isCurved){
                double rmin = m.stave_r + m.sensor_r;
                double half_width = abs(component.xmax()-component.xmin())/(2.*M_PI*rmin)*(2.0*M_PI)/2.;
                double phi_offset = getAttrOrDefault(component, _Unicode(phi_offset), double(0.0));
                Tube ele_box = Tube(rmin, rmin + m.sensor_thickness, abs(component.ymax()-component.ymin())/2., -half_width+phi_offset, half_width+phi_offset);
                ele_vol = Volume("sensor" +  _toString(iSensor, "_%d"), ele_box, m.sensor_material);                    
            }
            else{
                Box ele_box = Box(m.sensor_thickness/2., abs(component.xmax()-component.xmin())/2., abs(component.ymax()-component.ymin())/2.);            
                ele_vol = Volume("sensor" +  _toString(iSensor, "_%d"), ele_box, m.sensor_material);                    
            }

            if(m.sensor_sensitives.back())
                ele_vol.setSensitiveDetector(sens);
            ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());  
            m.sensor_volumes.push_back(ele_vol);

            iSensor++;
        }
        m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
        m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
        cout << det_name << ": Module: " << m.name << ", sensor width: " << to_string(m.sensor_width)  << ", sensor length: " << to_string(m.sensor_length) << endl;

        stave_information_list.push_back(m);
        cout << "Read stave information of stave " << m.name << endl;
    }

    int iModule_tot = 0;

    //=========  loop over layer elements in xml  ======================================

    cout << "Building " << det_name << " barrel detector " << det_name << "..." << endl;
    for(xml_coll_t c(e, _U(layer) ); c; ++c)  {

        xml_comp_t x_layer( c );

        // child elements: ladder, sensitive and periphery

        int layer_id = x_layer.id();
        int side = getAttrOrDefault(x_layer, _Unicode(side), 0);    // Use side=1 or -1 to use two staves/wafers with the same layer id

        double dr           = x_layer.dr(0);     // Spacing in r for every second stave. 
        double layer_offset = x_layer.offset(0);
        double z_offset     = x_layer.z_offset(0);

        string nameStr          = x_layer.nameStr();
        int    nmodules         = x_layer.nmodules();
        double step             = x_layer.step();   // Spacing of modules


        // Use the correct stave
        auto m = *find_if(stave_information_list.cbegin(), stave_information_list.cend(), [&nameStr] (const stave_information& stave) {
            return stave.name == nameStr;
        });    


        double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;
        
        double motherVolThickness = getAttrOrDefault(x_layer, _Unicode(motherVolThickness), double(5.0));
        double motherVolLength = getAttrOrDefault(x_layer, _Unicode(motherVolLength), double(stave_length));
        double motherVolOffset = getAttrOrDefault(x_layer, _Unicode(motherVolOffset), double(0.0)); // In case wafer/stave is asymmetric

        std::string layer_name = det_name+_toString(layer_id,"_layer%d")+_toString(side,"_side%d");
        double motherVolRmin = getAttrOrDefault(x_layer, _Unicode(motherVolRmin), double(x_layer.r()));
        Tube whole_layer_tube = Tube(motherVolRmin, motherVolRmin+motherVolThickness, motherVolLength/2.);

        Volume whole_layer_volume = Volume(layer_name, whole_layer_tube, theDetector.material("Air"));
        whole_layer_volume.setVisAttributes(theDetector, x_layer.visStr());
        pv = envelope.placeVolume(whole_layer_volume, Position(0., 0., z_offset));
        pv.addPhysVolID("layer", layer_id);

        DetElement layerDE( sdet , _toString(layer_id,"layer_%d")+_toString(side,"_side%d"), x_det.id() );
        layerDE.setPlacement( pv ) ;

        int nLadders = x_layer.attr<int>(  _Unicode(nLadders) ) ;
        double dphi = 2.*M_PI / double(nLadders);
        double phi0 = x_layer.phi0(0);

        //--------- loop over ladders ---------------------------
        for(int iStave=0; iStave<nLadders; ++iStave) {
            double layer_r = x_layer.r() + ((iStave%2 == 0) ? 0.0 : dr) + m.stave_r;  // Offset every second stave in r
            double phi = phi0 + iStave * dphi;
            RotationZYX rot( phi , 0, 0  ) ;

            double r = layer_r;
            double r_offset_component = layer_offset;

            double x_pos = r*cos(phi) - r_offset_component*sin(phi);
            double y_pos = r*sin(phi) + r_offset_component*cos(phi);
            double z_pos = motherVolOffset;
            Position stave_pos = Position(x_pos, y_pos, z_pos);
            
            string stave_name = layer_name + _toString(iStave,"_stave%d");

            PlacedVolume whole_stave_volume_placed;
            Volume whole_stave_volume_v;
            Assembly whole_stave_volume_a;
            if(m.motherVolThickness>0.0 && m.motherVolWidth>0.0){

                r = layer_r + m.motherVolThickness/2.;
                x_pos = r*cos(phi) - r_offset_component*sin(phi);
                y_pos = r*sin(phi) + r_offset_component*cos(phi);
                Box whole_stave_box = Box(m.motherVolThickness/2., m.motherVolWidth/2., stave_length/2.);
                whole_stave_volume_v = Volume(stave_name, whole_stave_box, theDetector.material("Air"));
                whole_stave_volume_v.setVisAttributes(theDetector, m.staveVis);
                whole_stave_volume_placed = whole_layer_volume.placeVolume(whole_stave_volume_v, Transform3D(rot,Position(x_pos, y_pos, z_pos)));
            }
            else{
                //// Use just assembly to avoid overlap between stave volume and other volumes
                whole_stave_volume_a = Assembly(stave_name);
                whole_stave_volume_placed = whole_layer_volume.placeVolume(whole_stave_volume_a, Transform3D(rot,stave_pos));
            }
            
            // Place components
            for(auto& component : m.components_vec){
                Assembly component_assembly(stave_name + "_" + component.name);
                if(m.motherVolThickness>0.0 && m.motherVolWidth>0.0)
                    pv = whole_stave_volume_v.placeVolume(component_assembly, Position(-m.motherVolThickness/2., 0., 0.));
                else
                    pv = whole_stave_volume_a.placeVolume(component_assembly);
                
                for(int i=0; i<int(component.thicknesses.size()); i++){
                    if(component.isCurved[i]){
                        double r_component_curved = 0.0; // component.r + component.rs[i]+ component.thicknesses[i]/2.;
                        r_offset_component = component.offset + component.offsets[i];
                        x_pos = r_component_curved*cos(phi) - r_offset_component*sin(phi);
                        y_pos = r_component_curved*sin(phi) + r_offset_component*cos(phi);
                        z_pos = component.z_offset + component.z_offsets[i] + motherVolOffset;
                        Position pos(x_pos, y_pos, z_pos);
                        component_assembly.placeVolume(component.volumes[i], Transform3D(rot,stave_pos).Inverse()*Transform3D(rot,pos));
                    }
                    else{
                        x_pos = component.r + component.rs[i] + component.thicknesses[i]/2.;
                        y_pos = component.offset + component.offsets[i];
                        z_pos = component.z_offset + component.z_offsets[i] + motherVolOffset;
                        Position pos(x_pos, y_pos, z_pos);
                        component_assembly.placeVolume(component.volumes[i], pos);
                    }
                }
                component_assembly->GetShape()->ComputeBBox(); 
            }

            // Place end of stave structures
            for(auto& endOfStave : m.endOfStaves_vec){
                Assembly endOfStave_assembly(stave_name + "_" + endOfStave.name);
                if(m.motherVolThickness>0.0 && m.motherVolWidth>0.0)
                    pv = whole_stave_volume_v.placeVolume(endOfStave_assembly, Position(-m.motherVolThickness/2., 0., 0.));
                else
                    pv = whole_stave_volume_a.placeVolume(endOfStave_assembly);

                for(int i=0; i<int(endOfStave.volumes.size()); i++){
                    std::vector<int> endOfStave_sides = {1,-1};
                    if(endOfStave.side[i] != 0)
                        endOfStave_sides = {endOfStave.side[i]};
                    for(auto& endOfStave_side : endOfStave_sides){  // Place it on both sides of the stave
                        if(endOfStave.isCurved[i]){
                            double r_component_curved = 0.0; // endOfStave.r + endOfStave.rs[i] + endOfStave.thicknesses[i]/2; // Correct for the fact that a tube element's origin is offset compared to the origin of a box
                            r_offset_component = endOfStave.offset + endOfStave.offsets[i];
                            x_pos = r_component_curved*cos(phi) - r_offset_component*sin(phi);
                            y_pos = r_component_curved*sin(phi) + r_offset_component*cos(phi);
                            z_pos = stave_length/2.+endOfStave.lengths[i]/2.+endOfStave.dzs[i] + endOfStave.z_offset + endOfStave.z_offsets[i];
                            Position pos = Position(x_pos, y_pos, z_pos*endOfStave_side+motherVolOffset);
                            endOfStave_assembly.placeVolume(endOfStave.volumes[i], Transform3D(rot,stave_pos).Inverse()*Transform3D(rot,pos));
                        }
                        else{
                            x_pos = endOfStave.r + endOfStave.rs[i] + endOfStave.thicknesses[i]/2.;
                            y_pos = endOfStave.offset + endOfStave.offsets[i];
                            z_pos = stave_length/2.+endOfStave.lengths[i]/2.+endOfStave.dzs[i] + endOfStave.z_offset + endOfStave.z_offsets[i]; 
                            Position pos(x_pos, y_pos, z_pos*endOfStave_side+motherVolOffset);
                            endOfStave_assembly.placeVolume(endOfStave.volumes[i], pos);
                        }
                    }
                }
                endOfStave_assembly->GetShape()->ComputeBBox(); 
            }

            // Place sensor
            for(int iModule=0; iModule<nmodules; iModule++){
                x_pos = m.sensor_r + m.sensor_thickness/2. + (iModule%2 == 0 ? 0.0 : m.stave_dr);
                y_pos = m.sensor_offset;
                z_pos = motherVolOffset + -(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step;
                Position pos(x_pos, y_pos, z_pos);
                    
                string module_name = stave_name + _toString(iModule,"_module%d");
                Assembly module_assembly(module_name);
                if(m.motherVolThickness>0.0 && m.motherVolWidth>0.0)
                    pv = whole_stave_volume_v.placeVolume(module_assembly, Position(-m.motherVolThickness/2., 0., 0.));
                else
                    pv = whole_stave_volume_a.placeVolume(module_assembly);
                pv.addPhysVolID("module", iModule_tot);

                DetElement moduleDE(layerDE,module_name,x_det.id());
                moduleDE.setPlacement(pv);


                // Place all sensor parts
                int iSensitive = 0;
                for(int i=0; i<int(m.sensor_volumes.size()); i++){
                    double r_component_curved = 0.0;
                    if(m.isCurved[i]){ // curved sensors
                        double phi_i = phi + 2.*(m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.)/(layer_r+m.sensor_r);
                        r_component_curved = m.sensor_thickness/2. + (iModule%2 == 0 ? 0.0 : m.stave_dr);
                        x_pos = r_component_curved*cos(phi_i) - r_offset_component*sin(phi_i);
                        y_pos = r_component_curved*sin(phi_i) + r_offset_component*cos(phi_i);
                        z_pos = motherVolOffset -(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step;
                        pos = Position(x_pos, y_pos, z_pos);

                        Position pos2(0., 0., m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.);
                        RotationZYX rot2(phi_i, 0, 0);
                        pv = module_assembly.placeVolume(m.sensor_volumes[i], Transform3D(rot,stave_pos).Inverse()*Transform3D(rot2, pos)*Translation3D(pos2));

                        if(m.sensor_sensitives[i]) { // Define as sensitive and add sensitive surface
                            string sensor_name = module_name + _toString(iSensitive,"_sensor%d");
                            pv.addPhysVolID("sensor", iSensitive);                
                            DetElement sensorDE(moduleDE,sensor_name,x_det.id());
                            sensorDE.setPlacement(pv);

                            Vector3D ocyl(+(r_component_curved+layer_r/2.), 0., 0.) ;  
                            dd4hep::rec::SurfaceType type = dd4hep::rec::SurfaceType::Sensitive;
                            type.setProperty(dd4hep::rec::SurfaceType::Cylinder,true);

                            double half_width = abs(m.sensor_xmax[i]-m.sensor_xmin[i])/(m.stave_r + m.sensor_r)/2.;
                            dd4hep::rec::VolCylinder surf(m.sensor_volumes[i], type, m.sensor_thickness/2., m.sensor_thickness/2., ocyl, half_width*2., +M_PI/2.-phi0+m.sensor_phi_offsets[i]);
                            volSurfaceList(sensorDE)->push_back(surf);
                            iSensitive++;
                        }
                    }
                    else{ // not curved
                        x_pos = 0.0;
                        y_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
                        z_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
                        Position pos2(x_pos, y_pos, z_pos);
                        pv = module_assembly.placeVolume(m.sensor_volumes[i], pos+pos2);

                        if(m.sensor_sensitives[i]) { // Define as sensitive and add sensitive surface
                            string sensor_name = module_name + _toString(iSensitive,"_sensor%d");
                            pv.addPhysVolID("sensor", iSensitive);                
                            DetElement sensorDE(moduleDE,sensor_name,x_det.id());
                            sensorDE.setPlacement(pv);

                            Vector3D u( 0. , 1. , 0. ) ;
                            Vector3D v( 0. , 0. , 1. ) ;
                            Vector3D n( 1. , 0. , 0. ) ;
                            VolPlane surf( m.sensor_volumes[i] , dd4hep::rec::SurfaceType::Sensitive , m.sensor_thickness/2. , m.sensor_thickness/2. , u,v,n );
                            volSurfaceList(sensorDE)->push_back(surf);
                            iSensitive++;
                        }
                    }
                }
                iModule_tot++;
                module_assembly->GetShape()->ComputeBBox();
            }
            if(m.motherVolThickness>0.0 && m.motherVolWidth>0.0){}
            else{
                whole_stave_volume_a->GetShape()->ComputeBBox();
            }
        }
    }

    pv.addPhysVolID( "system", x_det.id() ).addPhysVolID("side",0 )  ;

    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    return sdet;
}

DECLARE_DETELEMENT(VertexBarrel_o1_v02,create_element)
