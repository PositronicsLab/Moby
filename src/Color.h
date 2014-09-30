/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/


#include <osg/Array>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/NodeVisitor>
#include <osg/Vec4>
#include <osg/ShapeDrawable>
#include <osg/Texture2D>

namespace Moby {
/// This OSG NodeVisitor colors all Geometry & ShapeDrawable of a single node
 /**
  *  usage example:
  *    CcolorVisitor  newColor;
  *    newColor.setColor( r , g , b , a=1.0 );
  *    node->accept( newColor );
  **/
class  CcolorVisitor : public osg::NodeVisitor {

public :



    CcolorVisitor() : NodeVisitor( NodeVisitor::TRAVERSE_ALL_CHILDREN ) {
    // ---------------------------------------------------------------
    //
    // Default Ctors overide the default node visitor mode so all
    // children are visited
    //
    // ---------------------------------------------------------------

        //
        // Default to a white color
        //
        m_color.set( 1.0, 1.0, 1.0, 1.0 );

        m_colorArrays = new osg::Vec4Array;

        m_colorArrays->push_back( m_color );

        };


    CcolorVisitor( const osg::Vec4 &color ) : NodeVisitor( NodeVisitor::TRAVERSE_ALL_CHILDREN ) {
    // -------------------------------------------------------------------
    //
    // Overloaded Ctor initialised with the Color
    // Also override the visitor to traverse all the nodes children
    //
    // -------------------------------------------------------------------

        m_color = m_color;

        m_colorArrays = new osg::Vec4Array;

        m_colorArrays->push_back( m_color );


        };



    virtual ~CcolorVisitor(){};


    virtual
    void apply ( osg::Node &node ){
    // --------------------------------------------
    //
    //  Handle traversal of osg::Node node types
    //
    // --------------------------------------------

    traverse( node );

    } // apply( osg::Node &node )


    virtual
    void apply( osg::Geode &geode ){
    // ------------------------------------------------
    //
    //  Handle traversal of osg::Geode node types
    //
    // ------------------------------------------------
    osg::StateSet *state   = NULL;
    unsigned int    vertNum = 0;

    //
    //  We need to iterate through all the drawables check if
    //  the contain any geometry that we will need to process
    //
    unsigned int numGeoms = geode.getNumDrawables();


    for( unsigned int geodeIdx = 0; geodeIdx < numGeoms; geodeIdx++ ) {

        //
        // Use 'asGeometry' as its supposed to be faster than a dynamic_cast
        // every little saving counts
        //
       osg::Geometry *curGeom = geode.getDrawable( geodeIdx )->asGeometry();


        //
        // Only process if the drawable is geometry
        //
        if ( curGeom ) {
           osg::Vec4Array *colorArrays = dynamic_cast< osg::Vec4Array *>(curGeom->getColorArray());

           if ( colorArrays ) {
             for ( unsigned int i = 0; i < colorArrays->size(); i++ ) {

                osg::Vec4 *color = &colorArrays->operator [](i);

                //
                // could also use *color = m_color
                //
                color->normalize();
                color->set( color->_v[0]*m_color._v[0], color->_v[0]*m_color._v[1], color->_v[0]*m_color._v[2], color->_v[3]*m_color._v[3]);

             }
          } else{
              curGeom->setColorArray( m_colorArrays.get());
              curGeom->setColorBinding( osg::Geometry::BIND_OVERALL );
          }
      } else {
        // If Shape is not a geometry, it might be a shape, change shape color
        osg::ShapeDrawable * sd =  dynamic_cast<osg::ShapeDrawable*>(geode.getDrawable( geodeIdx ));
        if(sd){
          if(m_color._v[3] == 0.0){
            int s = 64;
            unsigned char* pData = new unsigned char[s*s*3];
            for ( int y=0; y<s; y++ )
              for ( int x=0; x<s; x++ )
              {
                unsigned char c = (( ((y&0x8)==0) ^ (((x&0x8))==0) ))*255;
                if(c > 127.5){
                  pData[x*3+y*s*3+0] = 255;
                  pData[x*3+y*s*3+1] = 255;
                  pData[x*3+y*s*3+2] = 255;
                }else{
                  pData[x*3+y*s*3+0] = 255*m_color._v[0];
                  pData[x*3+y*s*3+1] = 255*m_color._v[1];
                  pData[x*3+y*s*3+2] = 255*m_color._v[2];
                }
              }

            osg::Image* pImage = new osg::Image;
            pImage->setImage(
              s, s, 1,                    // 1=r? depth perhaps?
              GL_RGB,                     // internal format
              GL_RGB,GL_UNSIGNED_BYTE,		// pixelformat, type
              pData,                      // data
              osg::Image::USE_NEW_DELETE,	// mode
              1 );                        // packing

            osg::Texture2D* pTex = new osg::Texture2D;
            pTex->setImage( pImage );

            osg::StateSet* pStateSet = geode.getOrCreateStateSet();
            pStateSet->setTextureAttributeAndModes( 0, pTex, osg::StateAttribute::ON );
          } else {
            sd->setColor(osg::Vec4(m_color._v[0], m_color._v[1], m_color._v[2], m_color._v[3]));
          }
        }
      }
    }
  } // apply( osg::Geode


    void
    setColor( const float r, const float g, const float b, const float a = 1.0f ){
    // -------------------------------------------------------------------
    //
    // Set the color to change apply to the nodes geometry
    //
    // -------------------------------------------------------------------

        osg::Vec4 *c = &m_colorArrays->operator []( 0 );

        m_color.set( r,g,b,a );

        *c = m_color;


       } // setColor( r,g,b,a )


    void
    setColor( const osg::Vec4 &color  ){
    // -------------------------------------------------------------------
    //
    // Set the color to change apply to the nodes geometry
    //
    // -------------------------------------------------------------------

        osg::Vec4 *c = &m_colorArrays->operator []( 0 );

        m_color = color;

        *c = m_color;


       } // setColor( vec4 )





private :


    osg::Vec4 m_color;
    osg::ref_ptr< osg::Vec4Array > m_colorArrays;


 }; // class CcolorVisitor
}
