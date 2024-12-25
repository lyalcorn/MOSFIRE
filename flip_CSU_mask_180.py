import sys
import os

import random
import re
import xml.etree.ElementTree as ET

from astropy.io import fits
from astropy.table import Table, Column
from astropy import coordinates as c
from astropy import units as u
from astropy.time import Time

==#import mosfire
from mosfire.core import *
from mosfire.mask import *

def diff_from_center(leftnumber, rightnumber):
	left_distance=137.40-float(leftnumber)
	right_distance=137.40-float(rightnumber)
	return left_distance, right_distance
	
def find_newbarpos(left_distance, right_distance):
	leftbarpos=round(137.40-left_distance,3)
	rightbarpos=round(137.40-right_distance,3)
	return leftbarpos,rightbarpos
	
def opposing_number(num):
	if 1 <= int(num) <= 46:
		return 47 - int(num)
	else:
		raise ValueError("Input number must be between 1 and 46.")
		
def opposing_bar(num):
	if 1 <= int(num) <= 92:
		return 93 - int(num)
	else:
		raise ValueError("Input number must be between 1 and 92.")
	
def make_flippedguess(slitnumber, leftbarnumber, rightbarnumber, leftBarPositionMM, rightBarPositionMM, centerPositionArcsec, targetCenterDistance=''):
	centerPositionArcsec_flipped=-1*float(centerPositionArcsec)
	if targetCenterDistance!='':
		targetCenterDistance_flipped=-1*float(targetCenterDistance)
	leftbar_cen_distance, rightbar_cen_distance=diff_from_center(leftBarPositionMM,rightBarPositionMM)
	leftbar_cen_distance_flipped=-1*rightbar_cen_distance #left bar in flipped mask was originally right bar
	rightbar_cen_distance_flipped=-1*leftbar_cen_distance #ditto but with left and right switched
	leftBarPositionMM_flipped, rightBarPositionMM_flipped=find_newbarpos(leftbar_cen_distance_flipped, rightbar_cen_distance_flipped)
	slitnumber_flipped=opposing_number(slitnumber)
	leftbarnumber_flipped=opposing_bar(rightbarnumber)
	rightbarnumber_flipped=opposing_bar(leftbarnumber)
	if targetCenterDistance!='':
		flipped_dict={'slitNumber':str(slitnumber_flipped),'leftBarNumber':str(leftbarnumber_flipped), 'rightBarNumber':str(rightbarnumber_flipped),
					  'leftBarPositionMM':str(leftBarPositionMM_flipped), 'rightBarPositionMM':str(rightBarPositionMM_flipped), 
					  'centerPositionArcsec':str(centerPositionArcsec_flipped), 'targetCenterDistance':str(targetCenterDistance_flipped)}
		flipped_array=[slitnumber_flipped, leftbarnumber_flipped, rightbarnumber_flipped, leftBarPositionMM_flipped, rightBarPositionMM_flipped, centerPositionArcsec_flipped, targetCenterDistance_flipped]
	else:
		flipped_dict={'slitNumber':str(slitnumber_flipped),'leftBarNumber':str(leftbarnumber_flipped), 'rightBarNumber':str(rightbarnumber_flipped),
					  'leftBarPositionMM':str(leftBarPositionMM_flipped), 'rightBarPositionMM':str(rightBarPositionMM_flipped), 
					  'centerPositionArcsec':str(centerPositionArcsec_flipped)}
		flipped_array=[slitnumber_flipped, leftbarnumber_flipped, rightbarnumber_flipped, leftBarPositionMM_flipped, rightBarPositionMM_flipped, centerPositionArcsec_flipped]
	return flipped_dict
	
def flipmechslits(old_mask, root_xml):
	mechanicalSlitConfig = ET.SubElement(root_xml, 'mechanicalSlitConfig')
	transformed_mechslits_xml=[]
	for row1 in old_mask.slitpos:
		flipped_coordinates=make_flippedguess(row1['slitNumber'], row1['leftBarNumber'], 
											  row1['rightBarNumber'], row1['leftBarPositionMM'], 
											  row1['rightBarPositionMM'], row1['centerPositionArcsec'])
		transformed_mechslits_xml.append({'slitNumber':flipped_coordinates['slitNumber'],'leftBarNumber':flipped_coordinates['leftBarNumber'], 
										  'rightBarNumber':flipped_coordinates['rightBarNumber'], 'leftBarPositionMM':flipped_coordinates['leftBarPositionMM'], 
										  'rightBarPositionMM':flipped_coordinates['rightBarPositionMM'], 
										  'centerPositionArcsec':flipped_coordinates['centerPositionArcsec'], 'slitWidthArcsec':str(row1['slitWidthArcsec']), 
										  'target':str(row1['target'])})
	for slit_data in reversed(transformed_mechslits_xml):
		slit_element = ET.SubElement(mechanicalSlitConfig, 'mechanicalSlit')
		for key, value in slit_data.items():
			slit_element.set(key, value)
	return mechanicalSlitConfig
	
def flipscienceslits(old_mask, root_xml):
	scienceSlitConfig = ET.SubElement(root_xml, 'scienceSlitConfig')
	transformed_scienceslits_xml=[]
	for row1 in old_mask.scienceTargets:
		flipped_slitNumber=(len(old_mask.scienceTargets)+1)-int(row1['slitNumber'])
		flipped_targetCenterDistance=-1*float(row1['targetCenterDistance'])
		transformed_scienceslits_xml.append({'slitNumber':str(flipped_slitNumber),'slitRaH':str(row1['slitRaH']),
											 'slitRaM':str(row1['slitRaM']), 'slitRaS':str(row1['slitRaS']), 'slitDecD':str(row1['slitDecD']), 
											 'slitDecM':str(row1['slitDecM']), 'slitDecS':str(row1['slitDecS']), 
											 'slitWidthArcsec':str(row1['slitWidthArcsec']), 'slitLengthArcsec':str(row1['slitLengthArcsec']), 
											 'target':str(row1['target']), 'targetPriority':str(row1['targetPriority']), 
											 'targetMag':str(row1['targetMag']), 'targetCenterDistance':str(flipped_targetCenterDistance), 
											 'targetRaH':str(row1['targetRaH']), 'targetRaM':str(row1['targetRaM']), 
											 'targetRaS':str(row1['targetRaS']), 'targetDecD':str(row1['targetDecD']), 
											 'targetDecM':str(row1['targetDecM']), 'targetDecS':str(row1['targetDecS']), 
											 'targetEpoch':str(row1['targetEpoch']), 'targetEquinox':str(row1['targetEquinox'])})
	for slit_data in reversed(transformed_scienceslits_xml):
		slit_element = ET.SubElement(scienceSlitConfig, 'scienceSlit')
		for key, value in slit_data.items():
			slit_element.set(key, value)
	return scienceSlitConfig

def flipalignmentstars(old_mask, root_xml):
	alignment = ET.SubElement(root_xml, 'alignment')
	transformed_alignment_xml=[]
	for row1 in old_mask.alignmentStars:
		flipped_coordinates=make_flippedguess(row1['mechSlitNumber'], row1['leftBarNumber'], 
											  row1['rightBarNumber'], row1['leftBarPositionMM'], 
											  row1['rightBarPositionMM'], row1['centerPositionArcsec'], targetCenterDistance=row1['targetCenterDistance'])
		transformed_alignment_xml.append({'mechSlitNumber':flipped_coordinates['slitNumber'],'leftBarNumber':flipped_coordinates['leftBarNumber'], 
										  'rightBarNumber':flipped_coordinates['rightBarNumber'], 'leftBarPositionMM':flipped_coordinates['leftBarPositionMM'], 
										  'rightBarPositionMM':flipped_coordinates['rightBarPositionMM'], 
										  'centerPositionArcsec':flipped_coordinates['centerPositionArcsec'], 'slitWidthArcsec':str(row1['slitWidthArcsec']), 
										  'target':str(row1['target']), 'targetPriority':str(row1['targetPriority']), 'targetMag':str(row1['targetMag']), 
										  'targetCenterDistance':flipped_coordinates['targetCenterDistance'], 'targetRaH':str(row1['targetRaH']),
										  'targetRaM':str(row1['targetRaM']), 'targetRaS':str(row1['targetRaS']),
										  'targetDecD':str(row1['targetDecD']), 'targetDecM':str(row1['targetDecM']),
										  'targetDecS':str(row1['targetDecS']), 'targetEpoch':str(row1['targetEpoch']), 
										  'targetEquinox':str(row1['targetEquinox'])})
	for slit_data in reversed(transformed_alignment_xml):
		slit_element = ET.SubElement(alignment, 'alignSlit')
		for key, value in slit_data.items():
			slit_element.set(key, value)
	return alignment

def make_mascgenargs(old_mask, root_xml):
	mascgenElements=ET.SubElement(root_xml,'mascgenArguments')
	for tag, value in dict(old_mask.mascgenArguments).items():
		if isinstance(value, dict):
			# If the value is a dictionary, add it as attributes
			if tag=='maskName': 
				value=old_mask.name+'_180flip'
			tagname=ET.SubElement(mascgenElements, tag, attrib=value)
			for key, value2 in value.items():
				if key=='maskName': 
					value2=old_mask.name+'_180flip'
				tagname.set(key, value2)
			if tag=='inputs':
				centertags=ET.SubElement(tagname, 'center')
				centerattrib=old_mask.centerArguments
				for centerkey, centervalue in centerattrib['center'].items():
					centertags.set(centerkey, centervalue)
				steptags=ET.SubElement(tagname, 'steps')
				stepattrib=old_mask.stepArguments
				for stepkey, stepvalue in dict(stepattrib['steps']).items():
					steptags.set(stepkey, stepvalue)
			elif tag=='outputs':
				dirtags=ET.SubElement(tagname, 'directory')
				dirattrib=old_mask.dirArguments['directory']
				for dirkey, dirvalue in dict(dirattrib).items():
					dirtags.set(dirkey, dirvalue)
				scripttags=ET.SubElement(tagname, 'maskScript')
				scriptattrib=old_mask.scriptArguments['maskScript']
				for scriptkey, scriptvalue in dict(scriptattrib).items():
					scripttags.set(scriptkey, scriptvalue)
		else:
			# Otherwise, add it as text
			el = ET.SubElement(mascgenElements, tag)
			if tag=='maskName': 
				value=old_mask.name+'_180flip'
			el.text = str(value)
	return mascgenElements

def make_flipped_mask(filename, output_directory):
	print(f"Flipping CSU mask: {filename}")
	print(f"Output will be saved to: {output_directory}")
	
	old_mask=Mask(None)
	old_mask.read_xml(filename)
	output_filename_base=old_mask.name
	output_filename_xml=output_filename_base+'_180flip.xml'
	root = ET.Element('slitConfiguration')
	# Extract RA (Right Ascension) in hms
	ra_h, ra_m, ra_s = old_mask.center.ra.hms

	# Extract Dec (Declination) in dms
	dec_d, dec_m, dec_s = old_mask.center.dec.dms

	# Convert to strings
	ra_h_str = f"{int(ra_h):02}"  # RA Hour as zero-padded string
	ra_m_str = f"{int(ra_m):02}"  # RA Minute as zero-padded string
	ra_s_str = f"{ra_s:.2f}"      # RA Second as string with 2 decimal places

	dec_d_str = f"{int(dec_d):+03}"  # Dec Degree as signed string (e.g., +20, -01)
	dec_m_str = f"{int(dec_m):02}"   # Dec Minute as zero-padded string
	dec_s_str = f"{dec_s:.2f}"       # Dec Second as string with 2 decimal places
	maskDescription = ET.SubElement(root, 'maskDescription', maskName=old_mask.name+'_180flip', totalPriority=str(old_mask.priority), centerRaH=ra_h_str, centerRaM=ra_m_str, centerRaS=ra_s_str, centerDecD=dec_d_str, centerDecM=dec_m_str, centerDecS=dec_s_str, maskPA=str(float(old_mask.PA)+180.))

	mechanicalSlitConfig=flipmechslits(old_mask, root)

	scienceSlitConfig=flipscienceslits(old_mask, root)

	alignment=alignment=flipalignmentstars(old_mask, root)

	mascgenElements=make_mascgenargs(old_mask, root)
		
	tree = ET.ElementTree(root)
	tree.write(output_directory+output_filename_xml)
	
	print(f'Created:{output_directory}{output_filename_xml}')
	
if len(sys.argv) < 2:
	sys.exit(1)
filename = sys.argv[1]
output_directory = sys.argv[2] if len(sys.argv) > 2 else './'
if output_directory[-1]!='/': output_directory=output_directory+'/'

# Pass arguments to the main function
make_flipped_mask(filename, output_directory)

