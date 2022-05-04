#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces =
NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
    // Allocate memory for the quadric array
    size_t numVerts = mVerts.size();
    mQuadrics.reserve(numVerts);
    std::streamsize width = std::cerr.precision();  // store stream precision
    for (size_t i = 0; i < numVerts; i++) {

        // Compute quadric for vertex i here
        mQuadrics.push_back(createQuadricForVert(i));

        // Calculate initial error, should be numerically close to 0

        glm::vec3 v0 = mVerts[i].pos;
        glm::vec4 v(v0[0], v0[1], v0[2], 1);
        auto m = mQuadrics.back();

        // TODO CHECK
        auto error = glm::dot(v, (m * v));
        // std::cerr << std::scientific << std::setprecision(2) << error << " ";
    }
    std::cerr << std::setprecision(width) << std::fixed;  // reset stream precision

    // Run the initialize for the parent class to initialize the edge collapses
    DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute,
 * DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse* collapse) {
    // Compute collapse->position and collapse->cost here
    // based on the quadrics at the edge endpoints


    size_t v1 = e(collapse->halfEdge).vert;
    size_t v2 = e(e(collapse->halfEdge).pair).vert;

    //Q = Q1 + Q2, Q1 & Q2 = createQuadradicForVerts()
    glm::mat4 Q1 = createQuadricForVert(v1);
    glm::mat4 Q2 = createQuadricForVert(v2);
    glm::mat4 Q = Q1 + Q2;

    //To do this, we associate a symmetric 4×4 matrix Q with each vertex
    glm::mat4 Q_sym = Q;

    //The bottom row of the matrix is empty because v¯ is an homogeneous vector
    Q[0][3] = 0.f;
    Q[1][3] = 0.f;
    Q[2][3] = 0.f;
    Q[3][3] = 1.f;

    //If this also fails, we fall back on choosing v¯ from amongst the endpoints and the midpoint.
    //Check if Q is invertable
    const float EPSILON = 0.000015;

    if (glm::abs(glm::determinant(Q)) > EPSILON) {
        glm::vec4 v = glm::inverse(Q) * glm::vec4(0.f, 0.f, 0.f, 1.f);
        collapse->position = v;
        //glm does matrix multiplication in whatever way works, should not need to transpose v?
        collapse->cost = glm::dot(v , (Q_sym * v));
    } else {
        collapse->position = (mVerts[v1].pos + mVerts[v2].pos) / 2.f;
        collapse->cost = glm::length(collapse->position - mVerts[v1].pos);
    }
    
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(size_t ind) {
    DecimationMesh::updateVertexProperties(ind);
    mQuadrics[ind] = createQuadricForVert(ind);
}

/*! Lab 2
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
glm::mat4 QuadricDecimationMesh::createQuadricForVert(size_t indx) const {
    glm::mat4 Q({0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f});

    // The quadric for a vertex is the sum of all the quadrics for the adjacent
    // faces Tip: Matrix4x4 has an operator +=

    std::vector<size_t> neighborFaces = FindNeighborFaces(indx);

    for (int i = 0; i < neighborFaces.size(); i++) {
        Q += createQuadricForFace(neighborFaces.at(i));
    }

    return Q;
}

/*! Lab 2
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
glm::mat4 QuadricDecimationMesh::createQuadricForFace(size_t indx) const {

    // Calculate the quadric (outer product of plane parameters) for a face
    // here using the formula from Garland and Heckbert
    float a, b, c, d;

    //n^t = (a,b,c)
    glm::vec3 normal = glm::normalize(f(indx).normal);
    glm::vec3 point = v(e(f(indx).edge).vert).pos;

    a = normal[0];
    b = normal[1];
    c = normal[2];

    d = -1.f * glm::dot(point, normal); 


    glm::mat4 Kp = {{a*a, a*b, a*c, a*d}, 
                    {a*b, b*b, b*c, b*d}, 
                    {a*c, b*c, c*c, c*d},
                    {a*d, b*d, c*d, d*d} };

    return Kp;
}

void QuadricDecimationMesh::Render() {
    DecimationMesh::Render();

    glEnable(GL_LIGHTING);
    glMatrixMode(GL_MODELVIEW);

    if (mVisualizationMode == QuadricIsoSurfaces) {
        // Apply transform
        glPushMatrix();  // Push modelview matrix onto stack

        // Implement the quadric visualization here
        std::cout << "Quadric visualization not implemented" << std::endl;

        // Restore modelview matrix
        glPopMatrix();
    }
}
