macro(REGISTER_GLOBAL_C_SOURCES VARIABLE)
    file(GLOB C_SOURCES *.c)
    get_property(CURRENT GLOBAL PROPERTY ${VARIABLE})
    set_property(GLOBAL PROPERTY ${VARIABLE} "${CURRENT}" "${C_SOURCES}")
endmacro()
