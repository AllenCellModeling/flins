# API

What is the API for structures within the model?

## A base structure

Each site has several properties that are present:
 - x: the location within the tract that the site is located
 - parent: the parent protein or structure that this site is located on
 - binding site: the binding site that this structure uses to link to other sites
 - force: the force the site exerts if moved to an offset
 - energy: the energy necessary to move the site to an offset

 Each protein has several properties:
 - x: the location this protein is within a tract
 - length: how far the protein extends within the tract
 - force: the force the protein exerts if stretched 
