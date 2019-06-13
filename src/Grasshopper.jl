export GHType,
       GHInteger,
       GHNumber,
       GHPoint,
       GHBoolean

GHType(type, name, Nickname=name[1:1], Description=name, Access=:item, Default=0;
    nickname=Nickname, description=Description, access=Access, default=Default)=
    "$(type)|$(name)|$(nickname)|$(description)|$(access)|$(default)"
GHInteger(name, Nickname=name[1:1], Description=name, Access=:item, Default=0;
    nickname=Nickname, description=Description, access=Access, default=Default) =
    GHType("Integer", name, nickname, description, access, default)
GHNumber(name, Nickname=name[1:1], Description=name, Access=:item, Default=0;
    nickname=Nickname, description=Description, access=Access, default=Default) =
    GHType("Number", name, nickname, description, access, default)
GHPoint(name, Nickname=name[1:1], Description=name, Access=:item, Default=u0();
    nickname=Nickname, description=Description, access=Access, default=Default) =
    GHType("Point", name, nickname, description, access, default)
GHBoolean(name, Nickname=name[1:1], Description=name, Access=:item, Default=true;
    nickname=Nickname, description=Description, access=Access, default=Default) =
    GHType("Boolean", name, nickname, description, access, default)
