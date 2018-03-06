<?xml version='1.0'?>

<!-- XSL stylesheet configuration file for generating the documentation with Docbook-->
<!-- $Id: config.xsl,v 1.1 2006-09-08 16:02:44 nathan Exp $ -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
                version="1.0">

  <xsl:param name="use.id.as.filename" select="'1'"/>
  <xsl:param name="admon.graphics" select="'1'"/>
  <xsl:param name="admon.graphics.path"></xsl:param>
  <xsl:param name="chunk.section.depth" select="0"></xsl:param>

  <xsl:param name="passivetex.extensions" select="1"/>
  <xsl:param name="tex.math.in.alt" select="'latex'"/>
  
  <xsl:param name="chapter.autolabel" select="1"/>
  <xsl:param name="section.autolabel" select="1"/>

  <xsl:param name="html.stylesheet" select="'docbook.css'"/>

</xsl:stylesheet>