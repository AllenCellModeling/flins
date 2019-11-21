# encoding: utf-8
"""
Ged thinks such things are important. 

Proving a shorter unique name with which to reference objects in the world. 
"""

import uuid
import base64


def unique_name(name_uuid=None):
    """Create a unique name, from uuid if given"""
    if name_uuid is None:
        name_uuid = uuid.uuid4()
    name = base64.urlsafe_b64encode(name_uuid.bytes).strip(b"=").decode()
    return name


def name_to_uuid(name):
    """If, for some reason, you want the name back to UUID format"""
    name += "=" * (24 - len(name))  # repad stripped =
    name = name.encode()  # encode as bytes
    id = uuid.UUID(bytes=base64.urlsafe_b64decode(name))
    return id
