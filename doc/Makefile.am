if HAVE_DOXYGEN
    doxyfile.stamp:
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

    CLEANFILES = doxyfile.stamp
    all: doxyfile.stamp
    clean:
	rm -rf man latex html doxyfile.stamp
    distclean:
	rm -rf man latex html doxyfile.stamp Makefile
endif
