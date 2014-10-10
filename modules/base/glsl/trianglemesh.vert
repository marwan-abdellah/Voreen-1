/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2013 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

layout(location = 0) in vec4 vertexPosition;
#if defined(TRIANGLE_VEC4_VEC3)
layout(location = 1) in vec4 vertexColor;
layout(location = 2) in vec3 vertexNormal;
#elif defined(TRIANGLE_VEC3)
layout(location = 1) in vec3 vertexNormal;
#endif

out VertexData {
    vec4 color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
    vec3 normal;
#endif
} vertexData;

void main() {
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
    vertexData.normal = normalize(gl_NormalMatrix * vertexNormal);
#endif
#if defined(TRIANGLE_VEC4_VEC3)
    vertexData.color = vertexColor;
#else
    vertexData.color = vec4(1.0);
#endif
    gl_Position = vertexPosition;
}
